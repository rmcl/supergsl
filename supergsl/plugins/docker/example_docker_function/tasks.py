from invoke import task

@task
def hello(c):
    c.run('echo "hello world!"')

@task
def bash(c):
    c.run("bash", pty=True)
