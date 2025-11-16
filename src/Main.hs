module Main where

import Codec.Picture
import Control.Monad (forM_, unless)
import Control.Monad.ST
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as M

nx :: Int
nx = 64

ny :: Int
ny = 32

iters :: Int -- number of iterations
iters = 10

dt :: Double -- timestep
dt = 0.01

g :: Double -- gravity
g = -9.81

o :: Double -- overrelaxation
o = 1.9

h :: Double -- grid spacing
h = 1.0

density :: Double
density = 1.0

type MutVec = M.MVector RealWorld Double

data Sim = Sim
  { pressure :: !MutVec,
    smoke    :: !MutVec,
    newSmoke :: !MutVec,
    s        :: !MutVec,
    u        :: !MutVec,
    newU     :: !MutVec,
    v        :: !MutVec,
    newV     :: !MutVec
  }

data SampleType = U | V | S


idxS :: Int -> Int -> Int
idxS x y = x * ny + y

idxU :: Int -> Int -> Int
idxU x y = x * (ny + 1) + y

idxV :: Int -> Int -> Int
idxV x y = x * ny + y

coord :: Int -> (Int, Int)
coord i = (i `mod` nx, i `div` ny)

smokeColor :: Double -> Double -> PixelRGB8
smokeColor 0 _ = PixelRGB8 255 255 255
smokeColor _ smokeV =
  let c = floor (255 * smoke)
   in PixelRGB8 c c c

initSolids :: MutVec -> ST RealWorld ()
initSolids sV = forM_ [0 .. nx * ny - 1] $ \i -> do
  let (x, y) = coord i
  let isBorder = x == 0 || y == 0 || x == nx - 1 || y == nx - 1
  M.write sV i $ if isBorder then 0 else 1

diffuse :: MutVec -> MutVec -> ST RealWorld ()
diffuse vV sV = forM_ [0 .. nx * ny - 1] $ \i -> do
  oldV <- M.read vV i
  sCell <- M.read sV i
  if sCell == 0
    then return ()
    else M.write vV i $ oldV + g * dt

resetPressure :: MutVec -> ST RealWorld ()
resetPressure pressureV = forM_ [0 .. nx * ny - 1] $ \i -> do
  M.write pressureV i 0

project :: MutVec-> MutVec -> MutVec -> MutVec -> ST RealWorld ()
project uV vV sV pressureV = forM_ [0 .. iters] $ \_ -> do
  forM_ [0 .. nx * ny - 1] $ \i -> do
    let (x, y) = coord i
    sCell <- M.read sV i
    uCell <- M.read uV i
    vCell <- M.read vV i
    pCell <- M.read pressureV i

    unless (sCell == 0) $ do
      sx0 <- M.read sV $ idxS (x - 1) y
      sx1 <- M.read sV $ idxS (x + 1) y
      sy0 <- M.read sV $ idxS x (y - 1)
      sy1 <- M.read sV $ idxS x (y + 1)

      ux1 <- M.read uV $ idxU (x + 1) y
      vy1 <- M.read vV $ idxV x (y + 1)

      let d = ux1 - uCell + vy1 - vCell
          s' = sx0 + sx1 + sy0 + sy1
          pc = density * h / dt
          p = (-d) * o / s'
          p' = pCell + p * pc

      let u' = uCell - sx0 * p
          v' = vCell - sy0 * p
          u1' = ux1 + sx1 * p
          v1' = vy1 + sy1 * p

      M.write pressureV i p'
      M.write uV i u'
      M.write vV i v'
      M.write uV (idxU (x + 1) y) u1'
      M.write vV (idxV x (y + 1)) v1'

extrapolateBorder :: MutVec -> MutVec -> ST RealWorld ()
extrapolateBorder uV vV = do
  forM_ [0 .. nx - 1] $ \x -> do
    borderTopAdj <- M.read uV $ idxU x 1
    M.write uV (idxU x 0) borderTopAdj

    borderBottomAdj <- M.read uV $ idxU x (nx - 2)
    M.write uV (idxU x (nx - 1)) borderBottomAdj

  forM_ [0 .. ny - 1] $ \y -> do
    borderLeftAdj <- M.read vV $ idxV 1 y
    M.write vV (idxV 0 y) borderLeftAdj

    borderRightAdj <- M.read vV $ idxV (ny - 2) y
    M.write vV (idxV (ny - 1) y) borderRightAdj

copyVelocities :: MutVec -> MutVec -> MutVec -> MutVec -> ST RealWorld ()
copyVelocities newUV newVV uV vV = forM_ [1 .. nx * ny - 1] $ \i -> do
  u' <- M.read uV i
  v' <- M.read vV i
  M.write newUV i u'
  M.write newVV i v'

returnVelocities :: MutVec -> MutVec -> MutVec -> MutVec -> ST RealWorld ()
returnVelocities newUV newVV uV vV = forM_ [1 .. nx * ny - 1] $ \i -> do
  newU' <- M.read newUV i
  newV' <- M.read newVV i
  M.write uV i newU'
  M.write vV i newV'

advectVel :: MutVec -> MutVec -> MutVec -> ST RealWorld ()
advectVel uV vV sV = forM_ [0 .. nx * ny - 1] $ \i -> do
  let (x, y) = coord i
  unless (x == 0 || y == 0 || x == nx - 1 || y == ny - 1) $ do
    sCell <- M.read sV i
    sx0 <- M.read sV $ idxS (x - 1) y
    sy0 <- M.read sV $ idxS x (y - 1)
    uCell <- M.read uV i
    ux0 <- M.read uV $ idxU (x - 1) y
    ux1 <- M.read uV $ idxU (x + 1) y
    uy0 <- M.read uV $ idxU x (y - 1)
    uy1 <- M.read uV $ idxU x (y + 1)
    vCell <- M.read vV i
    vx0 <- M.read vV $ idxV (x - 1) y
    vx1 <- M.read vV $ idxV (x + 1) y
    vy0 <- M.read vV $ idxV x (y - 1)
    vy1 <- M.read vV $ idxV x (y + 1)

    unless (sCell == 0 || sx0 == 0) $ do
      -- check for the y > n cond
      let avgV = 0.25 * (vx0 + vx1 + vy0 + vy1)
      let xCoord = fromIntegral x * h - dt * uCell
          yCoord = fromIntegral y * h - dt * avgV
      nV <- sample vV xCoord yCoord V
      M.write vV i nV

    -- check for updating cell here
    unless (sCell == 0 || sy0 == 0) $ do
      let avgU = 0.25 * (ux0 + ux1 + uy0 + uy1)
      let xCoord = fromIntegral x * h - dt * avgU
          yCoord = fromIntegral y * h - dt * vCell
      nU <- sample uV xCoord yCoord U
      M.write uV i nU

sample :: MutVec -> Double -> Double -> SampleType -> ST RealWorld Double
sample vec xC yC field = do
  -- clamp to be within worldspace
  let xClamped = max 0.0 (min xC (fromIntegral nx - 1.0))
      yClamped = max 0.0 (min yC (fromIntegral ny - 1.0))

  -- scalars stored cell centres, velocities stored staggered
  let (offsetX, offsetY) =
        case field of
          U -> (0.0, 0.5 * h)
          V -> (0.5 * h, 0.0)
          S -> (0.5 * h, 0.5 * h)

  -- grid coord
  let x0 = floor (xClamped - offsetX / h) :: Int
      y0 = floor (yClamped - offsetY / h) :: Int

  let w00 = 1 - fromIntegral x0 / h
      w01 = fromIntegral x0 / h
      w10 = 1 - fromIntegral y0 / h
      w11 = fromIntegral y0 / h

  let idxX = case field of
        U -> idxU
        V -> idxV
        S -> idxS

  cell00 <- M.read vec $ idxX x0 y0
  cell01 <- M.read vec $ idxX x0 (y0 + 1)
  cell10 <- M.read vec $ idxX (x0 + 1) y0
  cell11 <- M.read vec $ idxX (x0 + 1) (y0 + 1)

  let val =
        w00 * w01 * cell00
          + w01 * w10 * cell10
          + w01 * w11 * cell01
          + w00 * w11 * cell11
  return val

copySmoke :: MutVec -> MutVec -> ST RealWorld ()
copySmoke smokeV newSmokeV = forM_ [1 .. nx * ny - 1] $ \i -> do
  s <- M.read smokeV i
  M.write newSmokeV i s

returnSmoke :: MutVec -> MutVec -> ST RealWorld ()
returnSmoke smokeV newSmokeV = forM_ [1 .. nx * ny - 1] $ \i -> do
  nS <- M.read newSmokeV i
  M.write smokeV i nS

advectSmoke :: MutVec -> MutVec -> MutVec -> MutVec -> ST RealWorld ()
advectSmoke uV vV sV newSmokeV = forM_ [0 .. nx * ny - 1] $ \i -> do
  let (x, y) = coord i
  sCell <- M.read sV i
  uCell <- M.read uV i
  vCell <- M.read vV i

  unless (sCell == 0) $ do
    ux1 <- M.read uV $ idxU (x + 1) y
    vy1 <- M.read vV $ idxV x (y + 1)

    let uC = uCell + ux1 * 0.5
        vC = vCell + vy1 * 0.5

    let xStep = fromIntegral x * h + 0.5 * h - dt * uC
        yStep = fromIntegral y * h + 0.5 * h - dt * vC

    newS <- sample newSmokeV xStep yStep S
    M.write newSmokeV i newS

simulate :: Int -> Sim -> ST RealWorld ()
simulate 0 _ = return ()
simulate i sim = do
  diffuse (v sim) (s sim)
  resetPressure (pressure sim)
  project (u sim) (v sim) (s sim) (pressure sim)
  extrapolateBorder (u sim) (v sim)

  -- copying and returning to prevent accessing updated values
  copyVelocities (newU sim) (newV sim) (u sim) (v sim)
  advectVel (u sim) (v sim) (s sim)
  returnVelocities (newU sim) (newV sim) (u sim) (v sim)

  copySmoke (smoke sim) (newSmoke sim)
  advectSmoke (u sim) (v sim) (s sim) (newSmoke sim)
  returnSmoke (smoke sim) (newSmoke sim)

  simulate (i - 1) sim

main :: IO ()
main = do
  simulation <- stToIO $ do
    uV         <- M.replicate ((nx + 1) * ny) 0 -- initialising a MAC grid
    vV         <- M.replicate (nx * (ny + 1)) 0
    newUV      <- M.replicate ((nx + 1) * ny) 0
    newVV      <- M.replicate (nx * (ny + 1)) 0
    sV         <- M.replicate (nx * ny) 0
    pressureV  <- M.replicate (nx * ny) 0
    smokeV     <- M.replicate (nx * ny) 0
    newSmokeV  <- M.replicate (nx * ny) 0
    sim <-
      Sim
        { u        = uV,
          v        = vV,
          newU     = newUV,
          newV     = newVV,
          s        = sV,
          pressure = pressureV,
          smoke    = smokeV,
          newSmoke = newSmokeV
        }

    initSolids (s sim)

    simulate 1 sim

    return sim

  frozen <- stToIO $ V.freeze simulation

  --let minP = V.minimum $ pressure frozen
  --    maxP = V.maximum $ pressure frozen

  let renderpx x y = smokeColor ((s frozen) V.! idxS x y) ((smoke frozen) V.! idxS x y)
  writePng "output.png" $ generateImage renderpx nx ny

-- REWRITE TODOS

-- helpers
-- main
-- init
-- integrate
-- project
-- advectvel
-- sample
-- advectsmoke
