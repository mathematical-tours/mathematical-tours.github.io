from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from common import *
from anisotropic_models import *
import struct
n=256;y,x=np.mgrid[:n,:n]/(n-1);f=simple_cartoon(x,y);p,tri=fitted_mesh(64);v=np.round(simple_cartoon(*p.T)*255).astype('uint8')
# Explicit connectivity preserves the anisotropic strip in decoding.
coords=np.round(p*65535).astype('<u2');cells=tri.astype('<u2');payload=b'P1M1'+struct.pack('<HHH',n,len(p),len(tri))+coords.tobytes()+v.tobytes()+cells.tobytes()
side,V,T=struct.unpack('<HHH',payload[4:10]);offset=10;q=np.frombuffer(payload,dtype='<u2',count=2*V,offset=offset).reshape(V,2)/65535;offset+=4*V
values=np.frombuffer(payload,dtype='uint8',count=V,offset=offset)/255;offset+=V;tri2=np.frombuffer(payload,dtype='<u2',count=3*T,offset=offset).reshape(T,3);g=interpolate_mesh(q,tri2,values,n)
jpeg,encoded,param=encode_image_at_budget(f,len(payload),'JPEG2000');budget=max(len(payload),len(encoded));out=ROOT/'figures'/Path(__file__).resolve().parent.relative_to(ROOT/'figures-code');out.mkdir(parents=True,exist_ok=True)
(out/'mesh.bin').write_bytes(payload+b'\0'*(budget-len(payload)));(out/'image.jp2').write_bytes(encoded);(out/'jpeg2000-stream.bin').write_bytes(encoded+b'\0'*(budget-len(encoded)))
fig,a=canvas(1,4,height=3.4,width=13);image(a[0,0],f,'Smooth piecewise cartoon');mesh_plot(a[0,1],q,tri2,f'{T} contour-aligned triangles');image(a[0,2],g,f'P1 decoding: {snr(f,g):.1f} dB');image(a[0,3],jpeg,f'JPEG2000: {snr(f,jpeg):.1f} dB');fig.suptitle(f'Equal transmitted budget: {8*budget/f.size:.3f} bits/pixel',fontsize=12)
finish(fig,__file__,'Uses smooth nonaffine intensities on both sides of an elliptical jump. A quasiuniform mesh in smooth regions transitions to thin, tangential triangles along the boundary; the same 1390-triangle mesh is used in Figures 5.20 and 5.21. The continuous P1 reconstruction is decoded from a real mesh stream containing a 10-byte header, uint16 coordinates, uint8 nodal values and uint16 triangle indices. The mesh and JPEG2000 streams have equal transmitted byte lengths, including all connectivity, headers and padding. No optimal entropy-coding claim is made.',parameters={'size':n,'angular_samples':64,'jpeg2000_parameter':param},checks={**mesh_checks(q,tri2),'vertices':V,'triangles':T,'equal_transmitted_bytes':budget,'mesh_payload_bytes':len(payload),'jpeg2000_payload_bytes':len(encoded),'P1_snr_db':snr(f,g),'jpeg2000_snr_db':snr(f,jpeg)},arrays={'points':q,'triangles':tri2,'P1_reconstruction':g,'jpeg2000_reconstruction':jpeg})
