function y = resample1(x,p,q)
y = uniformResample1(x,x,p,q);
end

function  [y, h] = uniformResample1(x, xTrue, p, q)

[p, q] = rat( p/q, 1e-12 );  %--- reduce to lowest terms
% (usually exact, sometimes not; loses at most 1 second every 10^12 seconds)
p = p(1);
q = q(1);
N = 10;
bta = 5;
pqmax = max(p,q);
fc = 1/2/pqmax;
L = 2*N(1)*pqmax + 1;
h = firls( L-1, [0 2*fc 2*fc 1], [1 1 0 0]).*kaiser(L,bta)' ;
h = p*h/sum(h);

Lhalf = (L-1)/2;
Lx = length(x);

nZeroBegin = floor(q-mod(Lhalf,q));
z = zeros(1,nZeroBegin);
h = [z h(:).'];  % ensure that h is a row vector.
Lhalf = Lhalf + nZeroBegin;

delay = floor(ceil(Lhalf)/q);
nZeroEnd = computeZeroPadLength(Lx,p,q,length(h),delay);
h = [h zeros(1,nZeroEnd)];

Ly = ceil(Lx*p/q);  % output length

szx = size(x);
sy = szx;
sy(1) = Ly;

sxTrue = size(xTrue);
syTrue = sxTrue;
syTrue(1) = Ly;

y = coder.nullcopy(complex(ones(syTrue),ones(syTrue)));

yVec = upfirdn(x,h,p,q);
indV = delay+(1:Ly);
indV = indV.';
yV = yVec(indV);
y = yV;
h = h(nZeroBegin+1:end-nZeroEnd);

end

function nZeroEnd = computeZeroPadLength(Lx,p,q,lenH,delay)
% coder.internal.prefer_const(Lx,p,q,lenH,delay);
nZeroEnd = 0;
while ceil( ((Lx-1)*p+lenH+nZeroEnd )/q ) - delay < ceil(Lx*p/q)
    nZeroEnd = nZeroEnd+1;
end
end