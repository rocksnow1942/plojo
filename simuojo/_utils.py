from os import path
import platform
systemname = platform.node()
systemname
if systemname.startswith('huis-mac-mini'):
    file_save_location = '/Users/hui/Cloudstation/R&D/Users/Hui Kang/Scripts/plojo/plojo/data'
else:
    file_save_location = path.join('/Users/hui/Scripts/plojo/plojo_app','data')
temp_position=path.join(file_save_location,'temp')
file_name = 'plojo_data'
