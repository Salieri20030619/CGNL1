import pathlib,importlib.util,subprocess,os,sys
R=pathlib.Path(r'C:\Users\31937\GC_validation_20260923');runtime=pathlib.Path(r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies');os.environ['PATH']=str(runtime/'native'/'poppler'/'Library'/'bin')+os.pathsep+os.environ['PATH']
skill=pathlib.Path(r'C:\Users\31937\.codex\plugins\cache\openai-primary-runtime\documents\26.909.12148\skills\documents\render_docx.py');spec=importlib.util.spec_from_file_location('skill_docx_renderer',skill);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
def convert(doc_path,user_profile,convert_tmp_dir,stem,verbose):
 out=pathlib.Path(convert_tmp_dir)/(stem+'.pdf');cmd=[r'C:\Windows\System32\WindowsPowerShell\v1.0\powershell.exe','-NoProfile','-ExecutionPolicy','Bypass','-File',str(R/'word_to_pdf.ps1'),'-DocPath',doc_path,'-PdfPath',str(out)]
 p=subprocess.run(cmd,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=180,creationflags=subprocess.CREATE_NO_WINDOW)
 debug=f'Native Word read-only PDF export, exit={p.returncode}\n{p.stdout}\n{p.stderr}'
 if p.returncode or not out.exists():return '',debug
 return str(out),debug
mod.convert_to_pdf=convert
sys.argv=['render_docx.py',str(R/'delivery'/'胃癌腹膜转移研究建议书.docx'),'--output_dir',str(R/'render_client_final'),'--dpi','144','--emit_pdf'];mod.main()
