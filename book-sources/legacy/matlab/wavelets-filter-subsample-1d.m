% Low/High pass filtering followed by sub-sampling.
a = subsampling( cconv(f,h) );
d = subsampling( cconv(f,g) );
% Up-sampling followed by filtering.
f1 =  cconv(upsampling(a),reverse(h)) + cconv(upsampling(d),reverse(g));
% Check that we really recover the same signal.
disp(strcat((['Error |f-f1|/|f| = ' num2str(norm(f-f1)/norm(f))])));