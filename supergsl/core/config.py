settings = {}
try:
    settings = json.loads(open('.supergsl.json', 'r'))
except Exception, error:
    print('Error loading config', error)
