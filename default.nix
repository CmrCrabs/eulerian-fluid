{ pkgs ? import <nixpkgs> {} }:
pkgs.mkShell {
  buildInputs = with pkgs.haskellPackages; [
    zlib
    ghc
    cabal-install
    imagemagick
    ffmpeg
  ];
}
