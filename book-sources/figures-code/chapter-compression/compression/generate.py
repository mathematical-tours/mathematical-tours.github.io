from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import *


import struct
f=flower(256);a,s=wavelet_array(f,'db2',4);T=.12;q=(np.sign(a)*np.floor(abs(a)/T)).astype('int16');decoded=np.sign(q)*(abs(q)+.5)*T;g=wavelet_inverse(decoded,s,'db2');idx=np.flatnonzero(q).astype('<u2');vals=q.ravel()[idx].astype('<i2')
payload=struct.pack('<4sIId',b'WVSP',a.size,len(idx),T)+idx.tobytes()+vals.tobytes();out=ROOT/'figures'/Path(__file__).resolve().parent.relative_to(ROOT/'figures-code');(out/'support-code.bin').write_bytes(payload)
fig,axs=canvas(1,4,height=3.4)
image(axs[0,0],f,'Original');image(axs[0,1],f[80:144,80:144],'Magnified crop');image(axs[0,2],(q!=0).astype(float),f'Support: {len(idx)} coefficients');wavelet_boundaries(axs[0,2],s);image(axs[0,3],g,f'Reconstruction: {snr(f,g):.1f} dB')
assert np.all((q!=0)==(abs(a)>=T));assert np.max(abs(a-decoded))<=T+1e-12
finish(fig,__file__,'Recomputes the orthogonal wavelet coefficients, dead-zone integer quantization, support, and midpoint reconstruction. The binary example stores a header, unsigned 16-bit indices, and signed 16-bit nonzero symbols; its measured length includes the header. It is an explicit support codec, not JPEG-2000. Zero is reconstructed as zero, and nonzero symbols as sign(q)(|q|+1/2)T.',data_sources=['data/flower.png'],parameters={'wavelet':'db2','levels':4,'threshold':T,'binary_header':'WVSP, N:uint32, M:uint32, T:float64'},checks={'support_size':len(idx),'encoded_bytes':len(payload),'bits_per_pixel':8*len(payload)/f.size,'snr_db':snr(f,g)},arrays={'quantized_coefficients':q,'reconstruction':g})
