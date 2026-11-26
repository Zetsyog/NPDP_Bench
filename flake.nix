{
  description = "NPDP Bench - non-serial polyadic dynamic programming benchmarks";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs =
    { self, nixpkgs }:
    let
      systems = [
        "x86_64-linux"
        "aarch64-linux"
      ];
      forAllSystems = f: nixpkgs.lib.genAttrs systems (system: f nixpkgs.legacyPackages.${system});
    in
    {
      devShells = forAllSystems (pkgs: {
        default = pkgs.mkShell {
          # stdenv already provides gcc (with OpenMP), make, coreutils, awk, sed
          packages = with pkgs; [
            bc # used by tools/finetune_pluto.sh
            clang-tools # clang-format, for make format / format-check
            valgrind
            gdb
          ];

          # The nix cc-wrapper silently drops -march=native by default, and
          # hardening flags (fortify, stack protector...) change the generated
          # code. Disable both so benchmark builds match a regular toolchain.
          NIX_ENFORCE_NO_NATIVE = 0;
          hardeningDisable = [ "all" ];

          # Pluto (polycc) is not packaged in nixpkgs (nixpkgs#pluto is an
          # unrelated Kubernetes tool), so it is taken from the environment.
          shellHook = ''
            if [ -z "''${PLUTO:-}" ]; then
              echo "NPDP Bench: PLUTO is not set, only the original (non-pluto) targets can be built."
              echo "            export PLUTO=/path/to/polycc to enable the *_pluto targets."
            fi
          '';
        };
      });

      # `nix flake check` runs the test suite (tests/run_tests.sh) and the
      # clang-format check; used by CI.
      checks = forAllSystems (pkgs: {
        default = pkgs.stdenv.mkDerivation {
          name = "npdp-bench-tests";
          src = self;
          nativeBuildInputs = [ pkgs.valgrind ];

          NIX_ENFORCE_NO_NATIVE = 0;
          hardeningDisable = [ "all" ];

          dontConfigure = true;
          buildPhase = ''
            runHook preBuild
            # LeakSanitizer needs ptrace, which the build sandbox forbids
            ASAN_OPTIONS=detect_leaks=0 make test
            runHook postBuild
          '';
          installPhase = "touch $out";
        };

        # Sources must match .clang-format (fix with `make format`).
        format = pkgs.runCommand "npdp-bench-format" { nativeBuildInputs = [ pkgs.clang-tools ]; } ''
          cd ${self}
          make format-check
          touch $out
        '';
      });
    };
}
