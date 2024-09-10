# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
#
# Astropy documentation build configuration file.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# Note that not all possible configuration values are present in this file.
#
# All configuration values have a default. Some values are defined in
# the global Astropy configuration which is loaded here before anything else.
# See astropy.sphinx.conf for which values are set there.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
# sys.path.insert(0, os.path.abspath('..'))
# IMPORTANT: the above commented section was generated by sphinx-quickstart, but
# is *NOT* appropriate for astropy or Astropy affiliated packages. It is left
# commented out with this explanation to make it clear why this should not be
# done. If the sys.path entry above is added, when the astropy.sphinx.conf
# import occurs, it will import the *source* version of astropy instead of the
# version installed (if invoked as "make html" or directly with sphinx), or the
# version in the build directory (if "python setup.py build_sphinx" is used).
# Thus, any C-extensions that are needed to build the documentation will *not*
# be accessible, and the documentation will not build correctly.

import datetime
import os
import sys
from configparser import ConfigParser
from importlib import import_module

try:
    from sphinx_astropy.conf.v1 import *  # noqa
except ImportError:
    print("ERROR: the documentation requires the sphinx-astropy package to be installed")
    sys.exit(1)

try:
    import tomllib
except ImportError:
    # Help users on older alphas
    import tomli as tomllib
from pathlib import Path

# Grab minversion from pyproject.toml
with (Path(__file__).parents[1] / "pyproject.toml").open("rb") as f:
    pyproject = tomllib.load(f)


# -- General configuration ----------------------------------------------------

# By default, highlight as Python 3.
# highlight_language = 'python3'

# If your documentation needs a minimal Sphinx version, state it here.
# needs_sphinx = '1.2'

# To perform a Sphinx version check that needs to be more specific than
# major.minor, call `check_sphinx_version("x.y.z")` here.
# check_sphinx_version("1.2.1")

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns += ["_templates", "notebooks/README.rst", "notebooks/Debug", "changes"]

# This is added to the end of RST files - a good place to put substitutions to
# be used globally.
rst_epilog += """
"""

# -- Project information ------------------------------------------------------

# This does not *have* to match the package name, but typically does
project = pyproject["project"]["name"]
author = ",".join(pyproject["project"]["authors"][0]["name"])
copyright = "{0}, {1}".format(datetime.datetime.now().year, pyproject["project"]["authors"])

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.

import_module(pyproject["project"]["name"])
package = sys.modules[pyproject["project"]["name"]]

# The short X.Y version.
version = package.__version__.split("-", 1)[0]
# The full version, including alpha/beta/rc tags.
release = package.__version__

# -- Options for HTML output --------------------------------------------------

# A NOTE ON HTML THEMES
# The global astropy configuration uses a custom theme, 'bootstrap-astropy',
# which is installed along with astropy. A different theme can be used or
# the options for this theme can be modified by overriding some of the
# variables set in the global configuration. The variables set in the
# global configuration are listed below, commented out.

# Add any paths that contain custom themes here, relative to this directory.
# To use a different custom theme, add the directory containing the theme.
# html_theme_path = ['themes]

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes. To override the custom theme, set this to the
# name of a builtin theme or the name of a custom theme in html_theme_path.
# html_theme = None
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]

html_theme_options = {
    "logotext1": "Sting",  # white,  semi-bold
    "logotext2": "ray",  # orange, light
    "logotext3": ":docs",  # white,  light
}

extensions += [
    "matplotlib.sphinxext.plot_directive",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "nbsphinx",
    "IPython.sphinxext.ipython_console_highlighting",
]

# Custom sidebar templates, maps document names to template names.
# html_sidebars = {}

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
# html_logo = ''

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = "images/stingray_logo.ico"

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
# html_last_updated_fmt = ''

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = "{0} v{1}".format(project, release)

# Output file base name for HTML help builder.
htmlhelp_basename = project + "doc"

# -- Options for LaTeX output -------------------------------------------------

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [("index", project + ".tex", project + " Documentation", author, "manual")]

# -- Options for manual page output -------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [("index", project.lower(), project + " Documentation", [author], 1)]

# Trust the links from doi.org, even if they might have Client errors or other minor issues
linkcheck_ignore = [r"https://doi.org/"]

# -- Options for the edit_on_github extension ---------------------------------

edit_on_github_branch = "main"


# -- Resolving issue number to links in changelog -----------------------------
github_issues_url = "https://github.com/{0}/issues/".format(
    pyproject["project"]["urls"]["repository"]
)

# -- Configuration for nbsphinx -----------------------------------------------
# disable notebook execution
nbsphinx_execute = "never"

# -- Generate DOI listing from Zenodo -----------------------------------------
import json
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass

ZENODO_API_ENDPOINT = "https://zenodo.org/api/records/"

# The “concept DOI” refers to all versions of Stingray.
# We'll use it to locate each of the specific versions in turn.
# See https://help.zenodo.org/#versioning for details.
CONCEPT_DOI = "10.5281/zenodo.1490116"


@dataclass
class Release(object):
    version: str
    doi: str

    @property
    def zenodo_url(self):
        return f"https://zenodo.org/record/{self.doi.split('.')[-1]}"

    @property
    def github_url(self):
        return f"https://github.com/StingraySoftware/stingray/releases/tag/{self.version}"

    @property
    def bibtex_url(self):
        return self.zenodo_url + "/export/hx"


params = urllib.parse.urlencode({"q": f'conceptdoi: "{CONCEPT_DOI}"', "all_versions": 1})
try:
    with urllib.request.urlopen(ZENODO_API_ENDPOINT + "?" + params) as url:
        data = json.loads(url.read().decode("utf-8"))
except urllib.error.URLError:
    data = {"hits": {"hits": []}}

releases = []
for rec in data["hits"]["hits"]:
    version = rec["metadata"]["version"]
    if version[0] != "v":
        continue
    doi = rec["metadata"]["doi"]
    releases.append(Release(version, doi))

with open("_zenodo.rst", "w") as f:
    if releases:
        f.write(".. list-table::\n")
        f.write("   :header-rows: 1\n\n")
        f.write("   * - Stingray Release\n")
        f.write("     - DOI\n")
        f.write("     - Citation\n")
        for r in sorted(releases, key=lambda r: r.version, reverse=True):
            f.write(f"   * - `{r.version} <{r.github_url}>`__\n")
            f.write(f"     - `{r.doi} <{r.zenodo_url}>`__\n")
            f.write(f"     - `[Link to BibTeX] <{r.bibtex_url}>`__\n")
    else:
        # The file needs to exist, but the text degrades gracefully if it is
        # empty.
        pass
