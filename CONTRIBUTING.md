# Contributing to fiasco

There are numerous ways to contribute to fiasco, including by
providing code and documentation, suggesting and discussing ideas,
submitting issues and bug reports, and engaging the broader scientific
community.

## Contributing code or documentation to fiasco

### Preliminaries

Before contributing to the fiasco code base, one must [**join
GitHub**](https://github.com/join?source=header-home).  A free account
will suffice for you to have unlimited public repositories.  If you
are new to [git](https://git-scm.com/), helpful resources include
documentation on [git
basics](https://git-scm.com/book/en/v2/Getting-Started-Git-Basics) and
an [interactive git
tutorial](https://try.github.io/levels/1/challenges/1).  You must also
[install
git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
locally on your computer.  We highly recommend getting familiar with
git by going through these tutorials or a [Software
Carpentry](https://software-carpentry.org/) workshop prior to making
code contributions.

### Forking and cloning fiasco

After creating your GitHub account, go to the [main
repository](https://github.com/wtbarnes/fiasco) and **fork a copy of
fiasco to your account**.

Next you must **clone your fork to your computer**.  Go to the
directory that will host your fiasco directory, and run one of the
following commands (after changing *your-username* to your username).
If you would like to use HTTPS (which is the default and easier to set
up), then run:

```ShellSession
git clone https://github.com/your-username/fiasco.git
```

SSH is a more secure option, but requires you to [set up an SSH
key](https://help.github.com/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent/) beforehand.  The equivalent SSH command is:

```ShellSession
git clone git@github.com:your-username/fiasco.git
```

After cloning, we must tell git where the development version of
fiasco is by running:

```ShellSession
git remote add upstream git://github.com/wtbarnes/fiasco.git
```

To check on which remotes exist, run `git remote -v`.  You should get
something like this:

```ShellSession
origin		git@github.com:namurphy/fiasco.git (fetch)
origin		git@github.com:namurphy/fiasco.git (push)
upstream	git@github.com:wtbarnes/fiasco.git (fetch)
upstream	git@github.com:wtbarnes/fiasco.git (push)
```

### Branches, commits, and pull requests

Before making any changes, it is prudent to update your local
repository with the most recent changes from the development
repository:

```ShellSession
git fetch upstream
```

Changes to fiasco should be made using branches.  It is usually best
to avoid making changes on your master branch so that it can be kept
consistent with the upstream repository.  Instead we can create a new
branch for the specific feature that you would like to work on:

```ShellSession
git branch *your-new-feature*
```

It is generally a good practice to choose a descriptive branch name.
Switch to your new branch by running:

```ShellSession
git checkout *your-new-feature*
```

After checking out your branch, let your fork of fiasco know about it
by running:

```ShellSession
git push --set-upstream origin *your-new-feature*
```

It is also useful to configure git so that only the branch you are
working on gets pushed to GitHub:

```ShellSession
git config --global push.default simple
```

Once you have set up your fork and created a branch, you are ready to
make edits to fiasco.

Go ahead and modify files with your favorite text editor.  Be sure to
include tests and documentation with any new functionality.  We also
recommend reading about [best practices for scientific
computing](https://doi.org/10.1371/journal.pbio.1001745).  fiasco uses
the [PEP 8 style guide for Python
code](https://www.python.org/dev/peps/pep-0008/) and the [numpydoc
format for
docstrings](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt)
to maintain consistency and readability.  New contributors should not
worry too much about precisely matching these styles when first
submitting a pull request, as further changes to the style can be
suggested during code review.

You may periodically commit changes to your branch by running

```ShellSession
git add filename.py
git commit -m "*brief description of changes*"
```

Committed changes may be pushed to the corresponding branch on your
GitHub fork of fiasco using

```ShellSession
git push origin *your-new-feature*
```

or, more simply,

```ShellSession
git push
```

Once you have completed your changes and pushed them to the branch on
GitHub, you are ready to make a pull request.  Go to your fork of
fiasco in GitHub.  Select "Compare and pull request".  Add a
descriptive title and some details about your changes.  Then select
"Create pull request".  Other contributors will then have a chance to
review the code and offer contructive suggestions.  You can continue
to edit the pull request by changing the corresponding branch on your
fiasco fork on GitHub.  After a pull request is merged into the code,
you may delete the branch you created for that pull request.
