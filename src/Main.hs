module Main where

import Codec.Picture
import Data.Vector.Mutable

n :: Int -- grid size
n = 256

iters :: Double -- number of iterations
iters = 10.0

dt :: Double -- timestep
dt = 0.01

g :: Double -- gravity
g = -9.81

o :: Double -- overrelaxation
o = 1.9

h :: Double -- grid spacing
h = 1.0

data Cell = Cell
  { density :: Double,
    u :: Double,
    v :: Double,
    u0 :: Double,
    v0 :: Double,
    s :: Int,
    pressure :: Double,
    solid :: Bool
  } deriving (Eq)

idx :: Int -> Int -> Int
idx x y = x + y * n

coord :: Int -> (Int, Int)
coord i = (i `mod` n, i `div` n)

color :: Cell -> PixelRGB8
color Cell {solid = True} = PixelRGB8 255 0 0
color _ = PixelRGB8 0 0 255

initSolids :: [Cell] -> [Cell]
initSolids = zipWith f [0 ..]
  where
    f i cell = cell {solid = isBorder, s = s'}
      where
        (x, y) = coord i
        isBorder = x == 0 || y == 0 || x == n - 1 || y == n - 1
        s' = if isBorder
            then 1
            else 0

diffuse :: [Cell] -> [Cell]
diffuse = map d
  where
    d cell
      | cell == cell {solid = True} = cell
      | otherwise = cell { v = v0 cell + g * dt }

--project :: [Cell] -> Int -> [Cell]
--project grid i
--  | i == 0 = grid
--  | otherwise = do
--
--      project grid (i-1)


divergence :: [Cell] -> Int -> Int -> Double
divergence grid x y = undefined

main :: IO ()
main = do
  let grid =
        replicate
          (n * n)
          Cell
            { density = 0,
              u = 0,
              v = 0,
              u0 = 0,
              v0 = 0,
              pressure = 0,
              solid = False
            }

  let grid' = initSolids grid

  --let grid  = diffuse grid

  let renderpx :: Int -> Int -> PixelRGB8
      renderpx x y = color $ grid' !! idx x y

  writePng "output.png" $ generateImage renderpx n n

-- DONE sim consts: N, dt, g, n (iterations), h (grid spacing), o (divergence)
-- DONE output png showing walls and density
-- DONE diffuse() (add gravity)
-- TODO project() (make incompressible)
-- TODO advect() (move velocity field)
-- TODO linear_solve()
