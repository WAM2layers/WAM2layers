# Developer's guide

This section explains how you can modify WAM2layers. For example, to use input
data from other sources. We encourage you to consider contributing your changes
so others can also benefit from your work. Therefore, the explanation of the
inner workings of WAM2layers are prefaced with instructions on how to
collaborate on GitHub.

```{admonition} Code of conduct
We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone. We pledge to act and
interact in ways that contribute to an open, welcoming, diverse, inclusive, and
healthy community.

For more information, see the
[contributor-covenant](https://www.contributor-covenant.org/version/2/1/code_of_conduct/).
```

## Collaborative workflow

Collaborative development of WAM2layers takes place on GitHub. You will need a
[GitHub account](https://github.com/join).

The process is roughly as follows:

1. You encounter a problem or question, or want to suggest an improvement or new feature.
1. You [Open(s) an
   issue](https://github.com/WAM2layers/WAM2layers/issues/new/choose) to discuss
   your idea.
1. Other developers comment on the issue and together you decide on the best way
   forward.
1. If the discussion results in a plan to change the code, you (or someone else)
   go ahead and make the necessary changes on a local copy of the code.
1. While you're working, you can create a draft pull requests (PR).
1. When you are happy with the changes, you mark your PR as "ready for review".
1. Every pull requests is reviewed by at least one core developer.
1. When all review comments are addressed and checks pass, the PR is merged.
1. At this point, the new code is part of the repository, but it is still "unreleased".
1. The core developers decide when to make a new release.
1. When the code is released, a new version incorporing your updates appears on
   Zenodo and PyPI.

### Forking the repository

Depending on your situation, you may have to fork the repository first. Forking
means you make a duplicate of the repository under your own GitHub account. The
advantage is that you will be able make changes to your fork without special
permissions. The disadvantage is that you have to keep your fork in sync with
the "upstream" repository, i.e. the original WAM2layers repo. Whenever the code
is updated there, your fork does not automatically get those updates too.

An alternative to forking a repository is to add you to the contributor team of
WAM2layers. This will give you the rights to push to the original repo. If you
feel this would make sense for your situation, talk to us.

You can make a fork by navigating to the [original
repo](https://github.com/WAM2layers/WAM2layers) and pressing the "Fork" button
(near the right top of the page).

### Cloning the repo (with good authentication)

Cloning a repo means obtaining a local copy. The version on GitHub is referred
to as the "remote" copy. Navigate to (your fork of) the WAM2layers repo and find
the big green button that says "code". Now, there are a few options.

Cloning with HTTPS works out of the box. However, whenever you want to "push"
changes back to the remote, you'll need to authenticate. This can be annoying.

Cloning with SSH is very convenient, but you have to [set up an SSH key
first](https://docs.github.com/en/authentication/connecting-to-github-with-ssh).
If you work on different machines, you can set multiple SSH keys. Despite the
initial efforts, in our experience this is the preferred way to authenticate
with GitHub.

There is also the option to use the GitHub command line interface. This is a
[separate software](https://cli.github.com/) provided by GitHub. We have little
experience with it.

Once you've decided how to authenticate, copy the corresponding url. Then open a
terminal and paste the url preceded with `git clone`. In case of using SSH:

```
# With fork
git clone git@github.com:Peter9192/WAM2layers.git

# Without fork
git clone git@github.com:WAM2layers/WAM2layers.git
```

### Development installation - additional tooling

You can install your local copy of the code by following the same steps as the
[regular installation](./installation.md). However, by adding `[develop]`, some
additional dependencies are installed that are convenient for developing.

```
pip install --editable .[develop]
```

The `--editable` flag means that any changes you make to the code will be
"effective immediately". When you run the model again, the updated code is used.

The `[develop]` option tells pip to install not just the code, but some
additional packages that are listed under the "develop" header in the file
`pyproject.toml`. These packages help with linting (checking your code against
syntax/style guides), automatic formatting, building documentation, running
tests, and publising the package on PyPI.

### Create a new branch

It is good practice to create a new "branch" for each feature you are adding.
For example, this will create (and switch to) a new branch called
"developer-documentation":

```
git checkout -b developer-documentation
```

You can use `git branch` to see all your branches, and
`git checkout <branchname>` (without the `-b`) to switch to another (existing)
branch. In case of the example, I would use this branch while working on this
documentation page. But if I was interrupting my work to fix a bug on another
part of the code, I would switch to another branch. This keeps all changes that
belong together, well, together.

### Make changes

Now you are ready to make the changes you need and test whether they work. Keep
in mind that your changes should not break the code for other people. If you're
adding new code, please also add relevant documentation, and add an update to
the changelog. If your introducing new dependencies, add them to
`pyproject.toml`. For more information on the inner workings of WAM2layers, see
below.

### Follow the style guide

Try to follow the Python style guide
([pep8](https://peps.python.org/pep-0008/)). This is where the developer tools
come in handy:

```
# flake8 can check your code against PEP8
flake8 .

# isort can automatically format import statements
isort .

# we use black for automatic formatting
black .
```

### Verify your changes

When you make changes to the code, you can easily verify the outputs by running the tests for the `preprocess`, `backtrack` and `visualize` workflows by executing the following command:

```py
pytest tests/test_workflow.py
```

If the outputs changed, you will receive an error message. It also provides you the instruction about how to update the reference files if you want to keep the new results.

### Commit and push your changes

Next, you can commit them and push your branch to the remote:

```
# see changed files or detailed diff
git status
git diff

# stage all updated files (-u) or individual files
git add -u
git add docs/develop.md

# commit the staged changes with a short description
git commit -m "updated develop documentation"

# push your branch to the remote the first time
git push -u origin developer-docs

# later commits can be pushed directly
git push
```

Frequently committing and pushing your code helps to keep a clear project
history and can save you from losing work. Also, others will be able to see what
you're working on, which can prevent double work and can be useful to ask for
confirmation or advise in an early stage, before you put in a lot of effort.

```{tip}
If you're using VScode, you can also use the [git
extension](https://code.visualstudio.com/docs/sourcecontrol/overview). It
provides a more visual experience. Similarly, you could also use [Github
desktop](https://desktop.github.com/).
```

#### Pre-commit

To help you make clean commits, we have configured this repository with pre-commit.
This is completely optional, but we recommend you to give it a try. You can run
pre-commit once by typing `pre-commit run` before committing some changes. If you
like this, and you want to do it automatically with each commit, you can type
`pre-commit install`.

### Open a (draft) pull request

When you first push your changes to the remote and navigate to the repository on
GitHub, you will see a popup that suggests to open a new pull request. If you
click this, you will see a comparison between your branch and the 'target'. If
you're working on a fork, you can set the target to the "upstream" repo, i.e.
the original WAM2layers repository. Usually the `main` branch will be the
target.

Add a clear title and description (these are often edited later), pointing to
the relevant issue(s) that the PR addresses. On the green button you can select
"Open draft PR".

When you open a pull request, some automatic checks are triggered. For example,
a preview of the documentation is build based on the code in the pull request.
You can look at the details of the checks to preview the documentation, or, if
they fail, what causes the issue(s).

### The review process

When you are happy with your changes, you can mark the PR "ready for review".
GitHub will automatically suggest some reviewers, or you can pick your own
choice. Anyway, every PR must be reviewed and approved by a core developer
before it is merged.

Reviewing others' PRs can be a good way to contribute to the community. GitHub
provides great functionality for commenting on individual lines of code, making
suggestions, et cetera. Please be welcoming and constructive in your feedback. A
good review involves verifying that the new code is clear, well documented, and
does what it's supposed to do, testing that the code (still) works, checking
that it conforms to the guidelines.

When the reviewer is happy with the code, they can approve the changes. At that
point, the code can be merged and becomes part of the "main" codebase. However,
a new release is required for it to become part of the PyPI and Zenodo
publications.

### Making a release

These instructions are intended for core developers.

You need an [account on PyPI](https://pypi.org/account/register/) with owner
rights on WAM2layers.

- Create a new branch, e.g. `prep_release_vXXX`
- Update "version" in pyproject.toml and in citation.cff (use [semantic
  versioning](https://semver.org/))
- Make sure the code is neatly formatted (see above)
- Make sure the changelog is up to date; update the "Unreleased" header with the
  new version, and add a new "Unreleased" header above
- Review and merge the release branch
- [Make a release on GitHub](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository)
- Make sure to pull the latest changes to main in your local copy of the repo.
- Publish to PyPI. In your terminal, in the base path of the repository, type:

  ```bash
  python -m build
  twine upload dist/*
  ```

- Check that the new release is successfully rendered on Zenodo and PyPI.
- Verify that you can install the new version with pip.

## The inner workings of WAM2layers

In this section we'll briefly go over some of the key insights required to
effectively modify and contribute code to WAM2layers.

### Structure of the repository

The structure of the repository follows that of a [standard Python
package](https://packaging.python.org/en/latest/tutorials/packaging-projects/).

Here is a very high-level overview of the role of the most important
files/folders:

```
.
├── pyproject.toml     --> Package metadata and dependencies. Used to `pip install` the package.
├── src/wam2layers/    --> This folder contains the source code of the model
│   ├── __init__.py    --> Needed for standard package structure
│   ├── analysis/      --> Code related to analyzing the input/output files
│   ├── cli.py         --> Code for the command line interface - entry point of WAM2layers
│   ├── preprocessing/ --> Code for preprocessing of raw input data
│   ├── tracking/      --> Actual moisture tracking code
│   └── utils/         --> Utilities that can be used in other parts of the code
├── docs               --> Sphinx project to create beautiful documentation
│   ├── _static        --> Place for e.g. images used in documentation
│   ├── conf.py        --> Configuration of the Sphinx project
│   └── index.md, etc. --> Actual documentation pages written in markdown
├── .readthedocs.yaml  --> Configuration for automatically publishing the documentation
├── scripts/           --> Example (download) scripts, not strictly part of the model.
├── tests/             --> Code to test the functions in src/wam2layers
├── misc/              --> Placeholder for random files such as logo
└── README.md          --> Shown on the GitHub landing page of WAM2layers
├── LICENSE            --> Enables others to re-use the code
├── CITATION.cff       --> Citation info; used by GitHub and Zenodo
├── CHANGELOG.md       --> Used to keep track of updates to the code
├── .gitignore         --> Tell git to ignore certain files - useful for keeping a clean repo
```

### GitHub integrations

The repository is set up to synchronize with Zenodo ([more
info](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content))
and ReadTheDocs ([more
info](https://docs.readthedocs.io/en/stable/connected-accounts.html)). Zenodo is
used to get a DOI to make the software citable. Every release will get its own
DOI, and there is also a "concept DOI" that applies to all versions. ReadTheDocs
automatically builds and hosts this documentation site.

### Click command line interface

The command line interface of WAM2layers is build with
[click](https://click.palletsprojects.com/en/8.1.x/).

In `pyproject.toml` the following configuration:

```toml
[project.scripts]
wam2layers = "wam2layers.cli:cli"
```

tells your computer that the command `wam2layers` should look for a function
called `cli` in the folder `wa2layers`, find the file `cli.py`. This function,
copied below, does nothing, but you will see that the docstring corresponds to
what is shown on the command line when you call `wam2layers`. The `@click.group`
decorator tells click that this function has subcommands, which are defined
below it. For example, the subcommand `backtrack` points to `backtrack_cli`,
which is an alias for the function called `cli` in
`wam2layers/tracking/backtrack.py` imported at the top of the file:

```py
import click
from wam2layers.tracking.backtrack import cli as backtrack_cli

@click.group()
def cli():
    """Command line interface to WAM2layers."""
    pass

cli.add_command(backtrack_cli, name="backtrack")
```

If you look at that function, you will see that it prints a welcome message and
then calls the function `run_experiment`. Note that it takes `config_file` an
input argument.

```{note}
It might be easier if all the methods related to the CLI are grouped together in
a single `cli.py` file. The reason why this is not done, is that by keeping a
`cli` function in each file, we can also run them as standalone scripts, e.g.
`python backtrack.py example-config.yaml`.
```

### Config file

The configuration file is written in [YAML](https://yaml.org/). YAML is a
human-friedly file format, very suitable for writting configuration files. A
YAML file follows the same structure as a (nested) dictorary in Python:
key-value pairs.

To parse the config, we read in the yaml file, but on top of that we use a
package called [Pydantic](https://docs.pydantic.dev/). Using Pydantic is very
similar to using Python dataclasses, for example:

```python
from dataclasses import dataclass

@dataclass
class Config:
    a: str
    b: int = 10

config = Config(a='test')
print(config)
# Config(a='test', b=10)

print(config.a)
# test
```

In this small example, we first define an object called `Config`. It has two
attributes: `a` and `b`. `a` should be a string, `b` an integer. `b` has a
default value of 10. When we call `config = Config(a='test')`, we are create an
instance of the class `Config` and assign it to a new variable called `config`.
Note that `b` will get the default value of `10`, as we didn't supply it.

So far, this is a standard Python dataclass. What Pydantic adds to this, is
validation. If you try to instantiate the config object above with
`Config(a=10)`, it will just work. With Pydantic, it will tell you that `a`
should be a string. The code changes to use Pydantic are minimal:

```
from pydantic import BaseModel

class Config(BaseModel):
    a: str
    b: int = 10
```

In `config.py`, all settings are documented with Python docstrings. They follow
the guidelines in [PEP 257](https://peps.python.org/pep-0257/). This
documentation is automatically imported in the Sphinx project (see below).

### Documentation

The documentation for WAM2layers uses
[Sphinx](https://www.sphinx-doc.org/en/master/), the de-facto standard for
documenting Python projects. By default, Sphinx documentation is written in
RestructuredText file format (.rst). However, WAM2layers uses
[MyST-NB](https://myst-nb.readthedocs.io/en/latest/), a Sphinx extension that
allows you to write documentation in the easier "markdown" format (.md). It can
also render Jupyter notebooks straight into the documentation. See
[here](https://jupyterbook.org/en/stable/reference/cheatsheet.html#links) for a
nice cheatsheet.

In addition, we use several plugins, configured in `docs/conf.py`. Importantly,
[autodoc](https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html) is
used to import docstrings from the main source code into the documentation.
[Napoleon](https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html#module-sphinx.ext.napoleon)
is used to have a more convenient format for the docstrings. [Sphinx ReadTheDocs
theme](https://sphinx-rtd-theme.readthedocs.io/en/stable/) is used for styling
of the documentation pages.

The documentation is hosted on [ReadTheDocs](https://readthedocs.org/). The
repository is setup to integrate with readthedocs, such that the documentation
automatically builds for all versions of the model, and also for pull requests.
So you can always preview the documentation on GitHub.

To build the documentation locally, from the base of the repository, type

```bash
sphinx-build -nW --keep-going -b html docs docs/_build/html
```

### Tests

The folder tests is intended for test code. Currently the content is very
minimal, but it would be good to add more tests in the future. Tests allow you
to check whether your code still works after you make some (perhaps seemingly
unrelated) changes. They can be set up to be run automatically when you make a pull request.

[Pytest](https://docs.pytest.org/en/7.2.x/) is used to run the tests. For example

```
pytest tests/
```

For more information, please have a look at the documentation of pytest.

### Typing

Since Python 3.5 type hints can be added to Python code.
These denote what the types of objects and function arguments should be.
A simple example:

```py
def add(a: int, b: int) -> float:
   return float(a+b)
```
The `add` function requires two arguments of the integer type, and will return a 
floating point number.
Adding types to code can help in detecting issues early, document what the function 
inputs should be, as well as make it easier for a developer's IDE to give useful 
suggestions.
More info on typing can be found [here](https://docs.python.org/3/library/typing.html).

To check if types are defined correctly, we use ["mypy"](https://mypy-lang.org/).
Mypy will go through all the source files, as well as imported packages, and check for
errors.
To run mypy do:

```sh
mypy src/ --install-types --ignore-missing-imports
```

`install-types` only needs to be run the first time, afterwards this is not needed.

If you are unsure about adding types, it is not required so you can leave them out.

### Python version support

In WAM2layers we follow Numpy's version support policy ([NEP-29](https://numpy.org/neps/nep-0029-deprecation_policy.html)).

This means that new releases will only support minor Python versions (e.g. 3.9, 3.10) which have been released within the last 42 months.
