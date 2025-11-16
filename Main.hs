module Main where

import Codec.Picture
import Control.Monad (forM_, unless)
import Control.Monad.ST
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as M

nx :: Int
ny :: Int
iters :: Int
dt :: Double
gravity :: Double
o :: Double
h :: Double
radius :: Double
gap :: Double
frames :: Int
density :: Double
obstacle :: Bool

-- ---------------------- PARAMETERS --------------------------
nx = 192 -- horizontal grid size
ny = 108 -- vertical grid size
iters = 50 -- number of gauss-seidel iterations, keep above ~10 for reasonably accurate results
dt = 1 / 220 -- smaller timestep -> more stable, keep as small as possible
gravity = -9.81 -- up is +ve, down is -ve, diffuse() is currently disabled for cooler visualisations
o = 1.9 -- overrelaxation, must be between 1 and 2, shouldnt need to be changed
h = 1.0 -- grid spacing, in "worldspace"
density = 1000.0 -- density of fluid, if you change it looks silly

obstacle = False -- whether or not to enable the obstacle (circle)
radius = 0.15 -- radius of obstacle if enabled
gap = 0.2 -- make sure gap >> radius

frames = 1 -- number of frames to render, and output to frames/frame_{i}.png
-- -------------------------------------------------------------

type MVec = M.MVector RealWorld Double
data Sim = Sim
  { pressure :: !MVec,
    smoke    :: !MVec,
    newSmoke :: !MVec,
    s        :: !MVec,
    u        :: !MVec,
    v        :: !MVec,
    newU     :: !MVec,
    newV     :: !MVec
  }

data SampleType = U | V | S

idx :: Int -> Int -> Int
idx x y = x + y * nx

coord :: Int -> (Int, Int)
coord i = (i `mod` nx, i `div` nx)

initSolids :: MVec -> ST RealWorld ()
initSolids sV = forM_ [0 .. nx * ny - 2] $ \i -> do
  let (x, y) = coord i
      isBorder = x == 0 || y == 0 || x == nx - 1 || y == ny - 1
      (xCr, yCr) = (gap * fromIntegral nx * h, 0.5 * fromIntegral ny * h)
      (xC, yC) = (fromIntegral x * h, fromIntegral y * h)

  let inCircle = dist xC yC xCr yCr <= (h * radius * fromIntegral ny)
      isObstacle = if obstacle then inCircle else False
  M.write sV i $ if isBorder || isObstacle then 0 else 1

dist :: Double -> Double -> Double -> Double -> Double
dist x1 y1 x2 y2 = sqrt ((x2-x1) * (x2-x1) + (y2-y1)*(y2-y1))

diffuse :: MVec -> MVec -> ST RealWorld ()
diffuse vV sV = forM_ [0 .. nx * ny - 2] $ \i -> do
  oldV <- M.read vV i
  sCell <- M.read sV i
  unless (sCell == 0) $ do
    M.write vV i $ oldV + gravity * dt

resetPressure :: MVec -> ST RealWorld ()
resetPressure pressureV = forM_ [0 .. nx * ny - 2] $ \i -> do
  M.write pressureV i 0

project :: MVec-> MVec -> MVec -> MVec -> ST RealWorld ()
project uV vV sV pressureV = forM_ [0 .. iters] $ \_ -> do
  forM_ [0 .. nx * ny - 2] $ \i -> do
    let (x, y) = coord i
    sCell <- M.read sV i
    uCell <- M.read uV i
    vCell <- M.read vV i
    pCell <- M.read pressureV i

    unless (sCell == 0) $ do
      sx0 <- M.read sV $ idx (x - 1) y
      sx1 <- M.read sV $ idx (x + 1) y
      sy0 <- M.read sV $ idx x (y - 1)
      sy1 <- M.read sV $ idx x (y + 1)

      let s' = sx0 + sx1 + sy0 + sy1
      unless (s' == 0.0) $ do
        ux1 <- M.read uV $ idx (x + 1) y
        vy1 <- M.read vV $ idx x (y + 1)

        let d = ux1 - uCell + vy1 - vCell
            pc = density * h / dt
            p = (-d) * o / s'
            p' = pCell + p * pc

        M.write pressureV i p'

        let u' = uCell - sx0 * p
            ux1' = ux1 + sx1 * p
            v' = vCell - sy0 * p
            vy1' = vy1 + sy1 * p

        M.write uV i u'
        M.write uV (idx (x + 1) y) ux1'
        M.write vV i v'
        M.write vV (idx x (y + 1)) vy1'

extrapolateBorder :: MVec -> MVec -> ST RealWorld ()
extrapolateBorder uV vV = do
  forM_ [0 .. nx - 1] $ \x -> do
    borderTopAdj <- M.read uV $ idx x 1
    M.write uV (idx x 0) borderTopAdj

    borderBottomAdj <- M.read uV $ idx x (ny - 2)
    M.write uV (idx x (ny - 1)) borderBottomAdj

  forM_ [0 .. ny - 1] $ \y -> do
    borderLeftAdj <- M.read vV $ idx 1 y
    M.write vV (idx 0 y) borderLeftAdj

    borderRightAdj <- M.read vV $ idx (nx - 2) y
    M.write vV (idx (nx - 1) y) borderRightAdj

copyVelocities :: MVec -> MVec -> MVec -> MVec -> ST RealWorld ()
copyVelocities newUV newVV uV vV = forM_ [0 .. nx * ny - 2] $ \i -> do
  u' <- M.read uV i
  v' <- M.read vV i
  M.write newUV i u'
  M.write newVV i v'

returnVelocities :: MVec -> MVec -> MVec -> MVec -> ST RealWorld ()
returnVelocities newUV newVV uV vV = forM_ [0 .. nx * ny - 2] $ \i -> do
  newU' <- M.read newUV i
  newV' <- M.read newVV i
  M.write uV i newU'
  M.write vV i newV'

advectVel :: MVec -> MVec -> MVec -> MVec -> MVec -> ST RealWorld ()
advectVel newUV newVV uV vV sV = forM_ [0 .. nx * ny - 2] $ \i -> do
  let (x, y) = coord i

  unless (x == 0 || y == 0 || x == nx - 1 || y == ny - 1) $ do

    sCell <- M.read sV i
    sx0 <- M.read sV $ idx (x - 1) y
    sy0 <- M.read sV $ idx x (y - 1)

    uCell <- M.read uV i
    ux1 <- M.read uV $ idx (x + 1) y
    uy0 <- M.read uV $ idx x (y - 1)
    ux1y0 <- M.read uV $ idx (x+1) (y-1)
    let avgU = 0.25 * (uCell + uy0 + ux1 + ux1y0)

    vCell <- M.read vV i
    vx0 <- M.read vV $ idx (x - 1) y
    vy1 <- M.read vV $ idx x (y + 1)
    vx0y1 <- M.read vV $ idx (x-1) (y+1)
    let avgV = 0.25 * (vCell + vx0 + vy1 + vx0y1)

    unless (sCell == 0 || sx0 == 0 || y > ny - 1) $ do
      let xC = fromIntegral x * h
          yC = fromIntegral y * h + (0.5 * h)
      nU <- sample uV (xC - dt * uCell) (yC - dt * avgV) U
      M.write newUV i nU

    unless (sCell == 0 || sy0 == 0 || x > nx - 1) $ do
      let xC = fromIntegral x * h + (0.5 * h)
          yC = fromIntegral y * h
      nV <- sample vV (xC - dt * avgU) (yC - dt * vCell) V
      M.write newVV i nV

sample :: MVec -> Double -> Double -> SampleType -> ST RealWorld Double
sample vec xC yC field = do
  let h1 = 1.0 / h

  let x = max h $ min xC (fromIntegral nx * h)
      y = max h $ min yC (fromIntegral ny * h)

  let (dx, dy) =
        case field of
          U -> (0.0, 0.5 * h)
          V -> (0.5 * h, 0.0)
          S -> (0.5 * h, 0.5 * h)

  let x0 = min (nx-1) $ floor ((x-dx)*h1)
      y0 = min (ny-1) $ floor ((y-dy)*h1)

  let tx = ((x-dx) - fromIntegral x0*h) *h1
      x1 = min (x0+1) (nx-1)
      ty = ((y-dy) -fromIntegral y0*h) * h1
      y1 = min (y0+1) (ny-1)
      sx = 1.0 - tx
      sy = 1.0 - ty

  let idxX = case field of
        U -> idx
        V -> idx
        S -> idx

  val00 <- M.read vec $ idxX x0 y0
  val10 <- M.read vec $ idxX x1 y0
  val11 <- M.read vec $ idxX x1 y1
  val01 <- M.read vec $ idxX x0 y1

  let val = sx * sy * val00 +
            tx * sy * val10 +
            tx * ty * val11 +
            sx * ty * val01
  return val

copySmoke :: MVec -> MVec -> ST RealWorld ()
copySmoke smokeV newSmokeV = forM_ [0 .. nx * ny - 2] $ \i -> do
  smk <- M.read smokeV i
  M.write newSmokeV i smk

returnSmoke :: MVec -> MVec -> ST RealWorld ()
returnSmoke smokeV newSmokeV = forM_ [0 .. nx * ny - 2] $ \i -> do
  nS <- M.read newSmokeV i
  M.write smokeV i nS

advectSmoke :: MVec -> MVec -> MVec -> MVec -> MVec -> ST RealWorld ()
advectSmoke uV vV sV smokeV newSmokeV = forM_ [0 .. nx * ny - 2] $ \i -> do
  let (x, y) = coord i
  sCell <- M.read sV i
  uCell <- M.read uV i
  vCell <- M.read vV i

  unless (sCell == 0) $ do
    ux1 <- M.read uV $ idx (x + 1) y
    vy1 <- M.read vV $ idx x (y + 1)

    let uC = uCell + ux1 * 0.5
        vC = vCell + vy1 * 0.5

    let xStep = fromIntegral x * h + 0.5 * h - dt * uC
        yStep = fromIntegral y * h + 0.5 * h - dt * vC

    newS <- sample smokeV xStep yStep S
    M.write newSmokeV i newS

windTunnel :: MVec -> MVec -> ST RealWorld ()
windTunnel uV smokeV = do
  let windU = 2.0
      smokeAmount = 1.0
      mid = ny `div` 2

  forM_ [mid-1 .. mid+1] $ \y -> do
    currentU <- M.read uV (idx 1 y)
    M.write uV (idx 1 y) (currentU + windU)
    currentSmoke <- M.read smokeV (idx 1 y)
    M.write smokeV (idx 1 y) (currentSmoke + smokeAmount)


main :: IO ()
main = do
  let len = nx * ny
  (uV, vV, newUV, newVV, sV, pressureV, smokeV, newSmokeV) <- stToIO $ do
    uV         <- M.replicate len 0
    vV         <- M.replicate len 0
    newUV      <- M.replicate len 0
    newVV      <- M.replicate len 0

    sV         <- M.replicate len 0
    pressureV  <- M.replicate len 0
    smokeV     <- M.replicate len 0
    newSmokeV  <- M.replicate len 0

    initSolids sV
    pure (uV, vV, newUV, newVV, sV, pressureV, smokeV, newSmokeV)

  forM_ [(0 :: Int)..frames] $ \i -> do
    stToIO $ do
      windTunnel uV smokeV
      --diffuse vV sV
      resetPressure pressureV
      project uV vV sV pressureV
      extrapolateBorder uV vV

      copyVelocities newUV newVV uV vV
      advectVel newUV newVV uV vV sV
      returnVelocities newUV newVV uV vV

      copySmoke smokeV newSmokeV
      advectSmoke uV vV sV smokeV newSmokeV
      returnSmoke smokeV newSmokeV


    frozenS     <- stToIO $ V.freeze sV
    frozenSmoke <- stToIO $ V.freeze smokeV
    frozenPressure <- stToIO $ V.freeze pressureV

    let minP = V.minimum frozenPressure
        maxP = V.maximum frozenPressure

    let renderpx ix iy =
          let x = ix
              y = (ny -1) - iy
          in sciSmokeColor (frozenS V.! idx x y) (frozenPressure V.! idx x y) minP maxP (frozenSmoke V.! idx x y)

    writePng ("frames/frame_" ++ show i ++ ".png") $ generateImage renderpx nx ny
    putStrLn ("rendered frame " ++ show i)


sciSmokeColor :: Double -> Double -> Double -> Double -> Double -> PixelRGB8
sciSmokeColor 0 _ _ _ _ = PixelRGB8 255 255 255
sciSmokeColor _ p minP maxP smokeV =
  let eps = 1e-4
      d = maxP - minP
      val = if d == 0 then 0.5 else (max minP (min p (maxP - eps)) - minP) / d
      m = 0.25
      num = floor (val / m) :: Int
      sc = (val - fromIntegral num * m) / m
      (r,g,b) = case num of
        0 -> (0, sc, 1)
        1 -> (0, 1, 1 - sc)
        2 -> (sc, 1, 0)
        _ -> (1, 1 - sc, 0)
      to8 x = floor (255 * x * min 1.0 smokeV)
  in PixelRGB8 (to8 r) (to8 g) (to8 b)
