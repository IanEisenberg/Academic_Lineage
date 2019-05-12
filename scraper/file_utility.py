import json

def save_json(data, file_path, pretty=True):
	with open(file_path, "w") as file:
		if pretty:
			json.dump(data, file, indent=4)
		else:
			json.dump(data, file)

def load_json(file_path):
	with open(file_path, "r") as file:
		return json.load(file)

def load_text(file_path):
	file = open(file_path, "r")
	text = file.read()
	file.close()
	return text