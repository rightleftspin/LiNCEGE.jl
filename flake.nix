{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.11";
    flake-parts.url = "github:hercules-ci/flake-parts";
    devshell.url = "github:numtide/devshell";
  };

  outputs = inputs @ {
    self,
    nixpkgs,
    flake-parts,
    ...
  }:
    flake-parts.lib.mkFlake {inherit inputs;} {
      systems = ["x86_64-linux"];
      imports = [inputs.devshell.flakeModule];

      perSystem = {
        config,
        pkgs,
        system,
        ...
      }: {
        devshells.default = {
          #env = [
          #  {
          #    name = "julia";
          #    value = "${pkgs.lib.getDev pkgs.openssl}";
          #  }
          #];
          packages = with pkgs; [
            julia-bin
            jupyter
          ];
        };
      };
    };
}
