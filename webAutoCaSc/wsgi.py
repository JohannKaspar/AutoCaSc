import sys
project_home = "/home/johann/AutoCaSc/AutoCaSc_project_folder/webAutoCaSc"

if project_home not in sys.path:
	sys.path.append(project_home)

from Dash_UI_styled import server as application

if __name__ == "__main__":
	application.run()
