from invoke import task

@task(default=True)
def sgsl(c):
    c.run('python supergsl/sgsl.py', pty=True)

@task
def typecheck(c):
    c.run('mypy --config-file mypy.ini supergsl')

@task(typecheck)
def test(c):
    c.run('nosetests supergsl --with-coverage --cover-package=supergsl --cover-xml', pty=True)

@task(typecheck)
def lint(c):
    c.run('pylint `ls -R|grep .py$|xargs`', pty=True)

@task
def bash(c):
    c.run("bash", pty=True)
