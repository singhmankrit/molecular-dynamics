{
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";

  outputs =
    {
      self,
      nixpkgs,
    }:
    let
      supportedSystems = [
        "x86_64-linux"
        "x86_64-darwin"
        "aarch64-linux"
        "aarch64-darwin"
      ];
      forAllSystems = nixpkgs.lib.genAttrs supportedSystems;
      pkgs = forAllSystems (system: nixpkgs.legacyPackages.${system});

    in
    {
      devShells = forAllSystems (
        system:
        let
          sps = pkgs.${system};
          pyflame = sps.python313Packages.buildPythonPackage rec {
            pname = "pyflame";
            version = "0.3.2";
            src = sps.fetchPypi {
              inherit pname version;
              hash = "sha256-j15RRngb3e84ezMXCyfPxb6Qf64BeVFttWSnI/MOUSE=";
            };
          };
        in
        {
          default = pkgs.${system}.mkShellNoCC {
            packages = with pkgs.${system}; [
              (python313.withPackages (ppkgs: [
                ppkgs.numpy
                ppkgs.matplotlib
                self.packages.${system}.alive
                pyflame
              ]))
              ffmpeg-headless # needed to make the animations
              mpv # for watching the generated videos
              flamegraph
            ];

            LD_LIBRARY_PATH = "${pkgs.${system}.libGL}/lib";
          };

        }
      );
      packages = forAllSystems (
        system:
        let
          sps = pkgs.${system};
        in
        {
          alive = sps.python313Packages.alive-progress.overrideAttrs ({
            postInstall = ''rm $out/LICENSE'';
          });
        }
      );
    };
}
