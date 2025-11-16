{ pkgs ? import <nixpkgs> {} }:
pkgs.mkShell {
  buildInputs = with pkgs; [
    haskellPackages.ghc
    haskellPackages.cabal-install
    zlib
    imagemagick
    ffmpeg
  ];
}
