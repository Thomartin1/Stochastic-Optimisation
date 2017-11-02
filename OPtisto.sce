
policy = [1.0 1.0 1.0 1.0 0.0 -1.0 -1.0 -1.0 -1.0; 1.0 1.0 1.0 1.0 0.0 -1.0 -1.0 -1.0 -1.0; 1.0 1.0 1.0 1.0 0.0 -1.0 -1.0 -1.0 -1.0; 1.0 1.0 1.0 1.0 0.0 -1.0 -1.0 -1.0 -1.0; 1.0 1.0 1.0 1.0 0.0 -1.0 -1.0 -1.0 -1.0; 1.0 1.0 1.0 1.0 0.0 -1.0 -1.0 -1.0 -1.0; 1.0 1.0 1.0 1.0 0.0 -1.0 -1.0 -1.0 -1.0; 1.0 1.0 1.0 1.0 0.0 -1.0 -1.0 -1.0 -1.0; 1.0 1.0 1.0 1.0 0.0 -1.0 -1.0 -1.0 -1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
function cost=simulation_mc(x0,policy,N)
  cost=0;
  for i=1:N do
    x=x0;
     for t=1:TF-1 do
      w=W(grand(1,1,'uin',1,2));
      x=x+policy(t,x)+w;
end
    cost=cost+(x-xref)^2
  end
  cost=cost/N;
endfunction

function cost=simulation_ex(x0,policy) // Exact computation with the law of W Wa=all_w(TF-1);
cost=0;
for i=1:size(Wa,'r') do
    x=x0;
    for t=1:TF-1 do
      x=x+policy(t,x)+Wa(i,t);
    end
    cost=cost+(x-xref)^2
  end
  cost=cost/size(Wa,'r');
endfunction

function W=all_w(n)
// generated all the possible (W_1,...,W_(TF-1)) if n==1 then
    W=[-1;1]
  else
    Wn=all_w(n-1);
    W=[-1*ones(size(Wn,'r'),1),Wn;1*ones(size(Wn,'r'),1),Wn];
  end
endfunction;

function costs=simulation_dp(policy)
// evaluation by dynamic programming with fixed policy Vs=ones(TF,length(X))*%inf;
// Bellman function at time TF
Vs(TF,:)=(X-xref) .^2;
// Compute final value functions
// Loop backward over time:
  for t=(TF-1):-1:1 do
    for x=1:10 do
      // loop on noises
      EV=0;
      for iw=1:size(W,"*") do
        next_state=x+policy(t,x)+W(iw);
        EV=EV+P(iw)*Vs(t+1,next_state);
      end
      Vs(t,x)=EV;
    end
end
  costs=Vs(1,:);
endfunction

print(simulation_dp(policy))