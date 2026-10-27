from pathlib import Path
O=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3')
p=O/'scripts/screen_liver.py';s=p.read_text(encoding='utf-8-sig').replace("c=a.CGNL1_log2CPM,s=2,cmap='viridis'","c=np.where(good,a.CGNL1_log2CPM,np.nan),s=2,cmap='viridis'").replace("set_title('CGNL1 log2(CPM+1)',fontsize=10)","set_title('CGNL1 log2(CPM+1), QC spots',fontsize=10)");p.write_text(s,encoding='utf-8')
p=O/'scripts/download_liver.py';s=p.read_text(encoding='utf-8').replace("print('ALL24 DOWNLOADS VERIFIED',flush=True)","print('VERIFIED '+str(len(jobs))+' DOWNLOADS',flush=True)");p.write_text(s,encoding='utf-8')
