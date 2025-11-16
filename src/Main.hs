module Main where

import Codec.Picture
import Control.Monad (forM_, unless)
import Control.Monad.ST
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as M

nx :: Int
nx = 64

ny :: Int
ny = 64

iters :: Int -- number of iterations
iters = 10

dt :: Double -- timestep
dt = 1 / 30

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
idxS x y = x + y * nx

idxU :: Int -> Int -> Int
idxU x y = x + y * (nx+1)

idxV :: Int -> Int -> Int
idxV x y = x + y * nx

coord :: Int -> (Int, Int)
coord i = (i `mod` nx, i `div` nx)


initSolids :: MutVec -> ST RealWorld ()
initSolids sV = forM_ [0 .. nx * ny - 1] $ \i -> do
  let (x, y) = coord i
  let isBorder = x == 0 || y == 0 || x == nx - 1 || y == ny - 1
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

    borderBottomAdj <- M.read uV $ idxU x (ny - 2)
    M.write uV (idxU x (ny - 1)) borderBottomAdj

  forM_ [0 .. ny - 1] $ \y -> do
    borderLeftAdj <- M.read vV $ idxV 1 y
    M.write vV (idxV 0 y) borderLeftAdj

    borderRightAdj <- M.read vV $ idxV (nx - 2) y
    M.write vV (idxV (nx - 1) y) borderRightAdj

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

advectVel :: MutVec -> MutVec -> MutVec -> MutVec -> MutVec -> ST RealWorld ()
advectVel newUV newVV uV vV sV = forM_ [0 .. nx * ny - 1] $ \i -> do
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
      M.write newVV i nV

    -- check for updating cell here
    unless (sCell == 0 || sy0 == 0) $ do
      let avgU = 0.25 * (ux0 + ux1 + uy0 + uy1)
      let xCoord = fromIntegral x * h - dt * avgU
          yCoord = fromIntegral y * h - dt * vCell
      nU <- sample uV xCoord yCoord U
      M.write newUV i nU

sample :: MutVec -> Double -> Double -> SampleType -> ST RealWorld Double
sample vec xC yC field = do
  -- clamp to be within worldspace
  let h1 = 1.0 / h

  let x = max h $ min xC (fromIntegral nx * h)
      y = max h $ min yC (fromIntegral ny * h)

  -- scalars stored cell centres, velocities stored staggered
  let (dx, dy) =
        case field of
          U -> (0.0, 0.5 * h)
          V -> (0.5 * h, 0.0)
          S -> (0.5 * h, 0.5 * h)

  -- grid coord
  let x0 = min (nx-1) $ floor ((x-dx)*h1)
      y0 = min (ny-1) $ floor ((y-dy)*h1)

  let tx = ((x-dx) - fromIntegral x0*h) *h1
      x1 = min (x0 + 1) (nx -1)
      ty = ((y-dy) -fromIntegral y0*h) * h1
      y1 = min (y0+1) (ny-1)
      sx = 1.0 - tx
      sy = 1.0 -ty

  let idxX = case field of
        U -> idxU
        V -> idxV
        S -> idxS

  val00 <- M.read vec $ idxX x0 y0
  val10 <- M.read vec $ idxX x1 y0
  val11 <- M.read vec $ idxX x1 y1
  val01 <- M.read vec $ idxX x0 y1

  let val = sx * sy * val00 +
            tx * sy * val10 +
            tx * ty * val11 +
            sx * ty * val01
  return val

copySmoke :: MutVec -> MutVec -> ST RealWorld ()
copySmoke smokeV newSmokeV = forM_ [1 .. nx * ny - 1] $ \i -> do
  smk <- M.read smokeV i
  M.write newSmokeV i smk

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
        yStep = fromIntegral y * h + 0.5 * h + dt * vC

    newS <- sample newSmokeV xStep yStep S
    M.write newSmokeV i newS

windTunnel :: MutVec -> MutVec -> ST RealWorld ()
windTunnel uV smokeV = do
  let windU = 0.5
      smokeAmount = 0.5

  forM_ [(ny `div` 3)..(ny - (ny `div` 3))] $ \y -> do
    currentU <- M.read uV (idxS 1 y)
    M.write uV (idxU 1 y) (currentU + windU)
    currentSmoke <- M.read smokeV (idxS 1 y)
    M.write smokeV (idxS 1 y) (currentSmoke + smokeAmount)

main :: IO ()
main = do
  (uV, vV, newUV, newVV, sV, pressureV, smokeV, newSmokeV) <- stToIO $ do
    uV         <- M.replicate ((nx + 1) * ny) 0
    vV         <- M.replicate (nx * (ny + 1)) 0
    newUV      <- M.replicate ((nx + 1) * ny) 0
    newVV      <- M.replicate (nx * (ny + 1)) 0
    sV         <- M.replicate (nx * ny) 0
    pressureV  <- M.replicate (nx * ny) 0
    smokeV     <- M.replicate (nx * ny) 0
    newSmokeV  <- M.replicate (nx * ny) 0

    initSolids sV
    pure (uV, vV, newUV, newVV, sV, pressureV, smokeV, newSmokeV)

  forM_ [0..30] $ \i -> do
    stToIO $ do
      windTunnel uV smokeV
      diffuse vV sV
      resetPressure pressureV
      project uV vV sV pressureV
      extrapolateBorder uV vV

      copyVelocities newUV newVV uV vV
      advectVel newUV newVV uV vV sV
      returnVelocities newUV newVV uV vV

      copySmoke smokeV newSmokeV
      advectSmoke uV vV sV newSmokeV
      returnSmoke smokeV newSmokeV

    frozenS     <- stToIO $ V.freeze sV
    frozenSmoke <- stToIO $ V.freeze smokeV
    frozenPressure <- stToIO $ V.freeze pressureV

    let minP = V.minimum $ frozenPressure
        maxP = V.maximum $ frozenPressure

    let renderpx x y =
          --smokeColor (frozenS V.! idxS x y) (frozenSmoke V.! idxS x y)
          sciColor (frozenS V.! idxS x y) (frozenPressure V.! idxS x y) minP maxP

    writePng ("frames/frame_" ++ show i ++ ".png") $
      generateImage renderpx nx ny
    putStrLn ("rendered frame " ++ show i)


-- ERRORS
-- when nx != ny checquering

sciColor :: Double -> Double -> Double -> Double -> PixelRGB8
sciColor 0 _ _ _ = PixelRGB8 255 0 0
sciColor _ p minP maxP =
  let eps = 1e-4
      d = maxP - minP
      v = if d == 0 then 0.5 else (max minP (min p (maxP - eps)) - minP) / d
      m = 0.25
      num = floor (v / m) :: Int
      s = (v - fromIntegral num * m) / m
      (r,g,b) = case num of
        0 -> (0, s, 1)
        1 -> (0, 1, 1 - s)
        2 -> (s, 1, 0)
        _ -> (1, 1 - s, 0)
      to8 x = floor (255 * x)
  in PixelRGB8 (to8 r) (to8 g) (to8 b)




smokeColor :: Double -> Double -> PixelRGB8
smokeColor 0 _ = PixelRGB8 255 0 0
smokeColor _ smokeV =
  let c = floor (255 * smokeV)
   in PixelRGB8 c c c
