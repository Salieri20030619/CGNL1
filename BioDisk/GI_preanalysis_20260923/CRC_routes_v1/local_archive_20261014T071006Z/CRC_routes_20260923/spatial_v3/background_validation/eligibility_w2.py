from common import *
D=Q/'ELIGIBILITY';D.mkdir(exist_ok=True);samples=json.loads((I/'SAMPLES.json').read_text());rows=[]
for s in samples:
 assert s['species']=='Homo sapiens';rows.append({'sample':s['sample'],'dataset':s['dataset'],'title':s['title'],'donor_label':s['patient_label_provisional'],'role':role(s),'species':s['species'],'colon_rectum':'not_individually_confirmed','treatment':'not_yet_verified','pairing':'C_and_L_not_paired; repeated_L8_collapsed','independence_from_spatial_cohort':'not_verified','input_complete':all((I/p).exists() for p in s['files'])})
x=pd.DataFrame(rows);assert x.input_complete.all();x.to_csv(D/'sample_eligibility.tsv',sep='\t',index=False);dump(D/'RESULT.json',{'status':'usable_for_exploratory_source_reference','samples':len(x),'clinical_inference_ready':False,'notes':['No inferential clinical comparisons until donor/treatment/site audit','GSE231559 overlap with MDACC spatial cohort is not excluded']})
