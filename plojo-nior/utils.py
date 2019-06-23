from os import path
import platform
systemname = platform.node()

file_name = 'hplc_data'

if systemname.startswith('huis-mac-mini'):
    _file_path = path.dirname(__file__)
    file_path = '/Users/hui/Cloudstation/R&D Backup/Plojo backup'# file path on my mac
    upload_data_folder = _file_path
    temp_position = _file_path
else:
    file_path = "C:\\Users\\aptitude\\Aptitude-Cloud\\R&D Backup\\Plojo backup"
    upload_data_folder = "C:\\Users\\aptitude\\Aptitude-Cloud\\R&D\\Shared Data\\Plojo\\!HPLC_UPLOAD"
    temp_position = "C:\\Users\\aptitude\\Aptitude-Cloud\\R&D Backup\\Plojo backup\\hplc_temp"
