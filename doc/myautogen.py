import cosmoslik.plugins, os
from textwrap import dedent

if not os.path.exists("plugins"): os.mkdir("plugins")

print "Found documentation for:"
for (name,cls,_) in cosmoslik.plugins.get_all_plugins():
    if cls.__doc__!=None: 
        print '  %s'%name
        with open('plugins/%s.rst'%'.'.join(name.split('.')[2:]),'w') as f: f.write(dedent(cls.__doc__))


