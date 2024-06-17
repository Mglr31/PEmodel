turbtest.NofRealisation=10;
plot_PE(rtest,["mu"])
%%
var10=exctract1pointSoundFluctuation(rturbinhomo,1.5,10);
var100=exctract1pointSoundFluctuation(rturbinhomo,1.5,100);
var200=exctract1pointSoundFluctuation(rturbinhomo,1.5,200);

%%
f=figure;
plot(var10.NormP2)
hold on
plot(var100.NormP2)
hold on
plot(var200.NormP2)
hold on
%%
figure
histogram(var10.NormA,20)
figure
histogram(var100.NormA,20)
figure
histogram(var200.NormA,20)
