import distutils.sysconfig as sysconfig
import importlib
import subprocess
import sys
import zipfile
from datetime import datetime
from pathlib import Path
from types import ModuleType
from typing import Generator, Iterable, Optional

import watermark

from hal.config import cfg


def gen_imports(globals_: dict) -> Generator:
    for name, val in globals_.items():
        if isinstance(val, ModuleType):
            base_module = val.__name__.split(".")[0]
            yield base_module
        else:
            try:
                yield val.__module__.split(".")[0]
            except AttributeError:
                pass


def zipdir(path: Path, zipf: zipfile.ZipFile, root: str = ""):
    for file_path in path.glob("**/*"):
        if "__pycache__" in file_path.parts:
            continue

        root = root or path.stem
        if file_path.is_file():
            relpath = file_path.relative_to(path)
            zipf.write(file_path, Path(root) / relpath)


def is_git_repo():
    try:
        # Run git command to check if we're in a git repository
        subprocess.check_output(
            ["git", "rev-parse", "--is-inside-work-tree"], stderr=subprocess.STDOUT
        )
        return True
    except subprocess.CalledProcessError:
        return False
    except FileNotFoundError:
        # Git command not found
        return False


def reproduce(
    globals_: dict,
    packages: Optional[Iterable[str]] = None,
    **watermark_kwargs,
) -> Path:
    script_path = Path(globals_["__file__"])
    output_path = script_path.parent / "output"
    output_path.mkdir(exist_ok=True, parents=True)

    combined = set(packages or []) | set(gen_imports(globals_))

    # Remove builtins
    combined -= {sys.builtin_module_names}

    # Remove standard library
    stdlib_pth = sysconfig.get_python_lib(standard_lib=True)
    stdlib = {p.stem.replace(".py", "") for p in Path(stdlib_pth).iterdir()}
    combined -= stdlib

    # Remove hal/ava
    combined -= {"hal", "builtins", "ava"}

    mark_kwargs = dict(
        author="{{ cookiecutter.author_name }}",
        current_time=True,
        current_date=True,
        timezone=True,
        updated=True,
        python=True,
        machine=True,
    )

    if is_git_repo():
        mark_kwargs.update(
            githash=True,
            gitbranch=True,
        )

    mark_kwargs.update(watermark_kwargs)
    mark = watermark.watermark(
        globals_=globals_,
        packages=",".join(sorted(combined)),
        **mark_kwargs,  # type: ignore
    )
    # get __version__ manually for selected packages
    versions = {}
    for package in packages or []:
        imported = importlib.import_module(package)
        try:
            versions[package] = imported.__version__
        except AttributeError:
            pass

    # write to mark string
    if versions:
        mark += "\n\n"
        mark += "Additional version information:\n"

    for package, version in versions.items():
        mark += f"\n{package}=={version}"

    # Run the command and capture the output
    freeze = subprocess.run(
        ["uv", "pip", "freeze", "--no-color"], capture_output=True, text=True
    )
    freeze_str = freeze.stdout

    # write to root freeze.txt file
    if freeze_str:
        freeze_file = cfg.root / "freeze.txt"
        with open(freeze_file, "w") as f:
            f.write("# uv pip freeze output generated at ")
            f.write(datetime.now().isoformat())
            f.write("\n")
            f.write(freeze_str)

            if freeze.stderr:
                f.write("\n# Error:\n")
                f.write(freeze.stderr)

    script_root = script_path.parent
    with zipfile.ZipFile(
        output_path / "_rpr.zip", "w", zipfile.ZIP_DEFLATED
    ) as rpr_zip:
        py_files = script_root.glob("**/*.py")
        for f in py_files:
            rpr_zip.write(f, f.relative_to(script_root))
        rpr_zip.writestr("watermark.txt", mark)

        rpr_zip.writestr("pip_freeze.txt", freeze_str)
        if freeze.stderr:
            rpr_zip.writestr("pip_freeze_error.txt", freeze.stderr)

        toolbox_dir = cfg.root / "ava" / "toolbox"
        if toolbox_dir.exists():
            zipdir(toolbox_dir, rpr_zip, root="_toolbox")

    return output_path
