from os import path
import platform

systemname = platform.node()

if systemname.startswith('huis-mac-mini'):
    file_save_location = '/Users/hui/Cloudstation/R&D/Users/Hui Kang/Scripts/plojo/plojo/data'
    temp_position=path.join(file_save_location,'temp')
else:
    file_save_location =  path.join(path.dirname(__file__),'data')
    temp_position = "C:\\Users\\aptitude\\Aptitude-Cloud\\R&D Backup\\Plojo backup\\temp"


file_name = 'plojo_data'
