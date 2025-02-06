import requests

url = "https://www.patricbrc.org/api/genome/100.11"

response = requests.get(url)

if response.status_code == 200:
	with open('patric_ids_metadata_100.11.json', 'w') as file:
		file.write(response.text)
	print("File saved as patric_ids_metadata_100.11.json")
else:
	print(f"Failed to download data. HTTP status code: {response.status_code}")
