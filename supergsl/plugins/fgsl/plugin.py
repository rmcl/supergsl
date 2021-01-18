from supergsl.plugins.docker import DockerPlugin


class fGSLPlugin(DockerPlugin):
    image_tag = 'fgsl'


if __name__ == '__main__':
    p = fGSLPlugin()
    p.build()
