from invoke import task

@task
def sgsl(c):
    c.run('python sgsl')

@task
def typecheck(c):
    c.run('mypy --config-file mypy.ini supergsl')

@task(typecheck)
def test(c):
    c.run('nosetests supergsl')

@task
def bash(c):
    c.run("bash", pty=True)
