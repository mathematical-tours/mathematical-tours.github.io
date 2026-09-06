from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
fig,axs=canvas(1,4,height=3.6,width=13);records=[];bpp=.5
for col,(name,f) in enumerate([('Flower',flower(256)),('Mandrill',mandrill(256))]):
    budget=int(f.size*bpp/8)
    for k,codec in enumerate(['JPEG','JPEG2000']):
        g,payload,quality=encode_image_at_budget(f,budget,codec);actual=8*len(payload)/f.size
        records.append({'image':name,'codec':codec,'target_bpp':bpp,'bytes':len(payload),'actual_bpp':actual,'parameter':quality,'snr_db':snr(f,g)})
        image(axs[0,2*col+k],g,f'{name}: {codec}\n{actual:.3f} bpp; {snr(f,g):.1f} dB')
        out=ROOT/'figures'/Path(__file__).resolve().parent.relative_to(ROOT/'figures-code');out.mkdir(parents=True,exist_ok=True)
        (out/(name.lower()+'-'+codec+('.jpg' if k==0 else '.jp2'))).write_bytes(payload)
finish(fig,__file__,'One row compares JPEG and JPEG-2000 on the flower and a different second image, the author-requested mandrill. Both encoders share the same uncompressed reference and a 0.5 bits/pixel target. Labels count actual stream bytes including headers; panels are decoded from the saved streams.',parameters={'target_bits_per_pixel':bpp,'encoder':'Pillow with libjpeg/OpenJPEG'},checks={'encoding_results':records})
