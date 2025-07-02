This is a project to develop a photometry pipeline for JWST NIRCam data processing.

You can clone the project using the following command:
```
git clone https://github.com/Nikhil0504/Photometry_Pipeline.git
```

and install the required dependencies by running the following command:
```
uv sync
```
(make sure to have uv installed)

You can bump the version using the following commands:
```sh
$ uvx --from just-bin just
Available recipes:
    bump-and-tag type # Internal recipe to bump version and create git tag
    bump-major        # Bump version by major and create git tag
    bump-minor        # Bump version by minor and create git tag
    bump-patch        # Bump version by patch and create git tag
    default           # List available commands
    push-all          # Push both commits and tag to remote
    push-tag          # Push the latest version tag to remote
    tag-version       # Create git tag from current version if it doesn't exist
    version           # Show current version
```

Scripts:
You can access the webbpsf's script to make webbpsfs using:
```bash
$ uv run webbpsfs file.fits F200W output_folder psf_out

$ uv run webbpsfs -h # for getting the other optional parameters
```
