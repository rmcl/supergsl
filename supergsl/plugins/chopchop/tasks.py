from invoke import task

@task
def hello(c):
    c.run('echo "helloooo world!"', pty=True)

@task
def build_genome(c, input_genome_fasta_path):
    """Create a bowtie genome file for use by chopchop"""
    c.run('./bowtie/bowtie-build cenpk1137d.fa cenpk1137d')

@task
def bash(c):
    c.run("bash", pty=True)
