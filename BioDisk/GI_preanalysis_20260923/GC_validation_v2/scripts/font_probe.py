import pathlib,subprocess,json,importlib.util
print(subprocess.run(['fc-match','-f','%{file}',':lang=zh'],capture_output=True,text=True).stdout)
print('document_modules',[(n,importlib.util.find_spec(n) is not None) for n in ['docx','markdown','pdf2image']])
