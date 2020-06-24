tFuns = [];

tFuns(1).name="exponential";
tFuns(1).lt=@(s)1/(1+s);
tFuns(1).ilt=@(x)exp(-x);
tFuns(1).xvals=0.1:0.01:5;
tFuns(1).yrange=[-0.5,1.5];

tFuns(2).name="sine";
tFuns(2).lt=@(s)1/(1+s^2);
tFuns(2).ilt=@(x)sin(x);
tFuns(2).xvals=0.1:0.05:15;
tFuns(2).yrange=[-1.2,1.2];

tFuns(3).name="heavyside";
tFuns(3).lt=@(s)exp(-s)/s;
tFuns(3).ilt=@(x)floor(x>1);
tFuns(3).xvals=0.1:0.01:3;
tFuns(3).yrange=[-0.5,1.5];

tFuns(4).name="expheavyside";
tFuns(4).lt=@(s)exp(-s)/(1+s);
tFuns(4).ilt=@(x)floor(x>1)*exp(1-x);
tFuns(4).xvals=0.1:0.01:5;
tFuns(4).yrange=[-0.5,1.2];

tFuns(5).name="squarewave";
tFuns(5).lt=@(s)(1/s)/(exp(s)+1);
tFuns(5).ilt=@(x)mod(floor(x),2);
tFuns(5).xvals=0.1:0.02:10;
tFuns(5).yrange=[-1,2];

tFuns(6).name="staircase";
tFuns(6).lt=@(s)(1/s)/(exp(s)-1);
tFuns(6).ilt=@(x)floor(x);
tFuns(6).xvals=0.1:0.01:5;
tFuns(6).yrange=[-1,5];


% Set the parameters of the DEMO here:
testFunction = "squarewave";
funEvals = 100;

% compute the inverse and plot the comparison
ix = 1;
while(tFuns(ix).name~=testFunction)
    ix = ix + 1;
end

cme = matlab_ilt(tFuns(ix).lt, tFuns(ix).xvals, funEvals, "cme");
euler = matlab_ilt(tFuns(ix).lt, tFuns(ix).xvals, funEvals, "euler");
gaver = matlab_ilt(tFuns(ix).lt, tFuns(ix).xvals, funEvals, "gaver");
exact = arrayfun(tFuns(ix).ilt, tFuns(ix).xvals)';

plot([exact, cme, euler, gaver])
ylim(tFuns(ix).yrange)
legend("exact", "cme", "euler", "gaver")

