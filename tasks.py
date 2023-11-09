from invoke import task

@task(default=True)
def sgsl(c):
    """Start the sgsl shell."""
    c.run('python supergsl/sgsl.py', pty=True)


@task
def typecheck(c):
    """Run the mypy static analysis."""
    c.run('mypy --config-file mypy.ini supergsl')


@task(typecheck)
def test(c):
    "Run nostests and generate coverage report."""
    c.run('PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python python -m pytest supergsl', pty=True)


@task(typecheck)
def lint(c):
    """Run the linter."""
    c.run('pylint `ls -R|grep .py$|xargs`', pty=True)


@task
def bash(c):
    """Invoke a bash shell."""
    c.run("bash", pty=True)


@task(test)
def build_dist(c):
    """Steps to build and package supergsl as library for pypi."""
    c.run('pip install twine')
    c.run('rm -rf build dist')
    c.run('python setup.py sdist bdist_wheel')
    c.run('twine check dist/*')
    c.run('echo "Build complete, run \'twine upload dist/*\' to upload to pypi."')
