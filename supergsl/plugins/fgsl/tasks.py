from invoke import task

@task
def gsl(c):
    c.run('mono Gslc/bin/Gslc/Gslc.exe', pty=True)

@task
def bash(c):
    c.run("bash", pty=True)
