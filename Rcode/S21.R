input = read.csv('data/S21_5.csv')
f = input$f*1e9
s21 = input$s21
h0 = 10^(input$s21[1]/(20))
Hf2 = (10^(input$s21/(20))/h0)^2
data = data.frame(f = f, Hf2 = Hf2)

P0 = 1.904172485e-3
q = 1.6e-19
ita = 3.645514e-01
I = 7.5e-3
Ith = I-P0/ita
# y = 186517947828
# z = 1.514458e22

#### functions #####
Ns_pre = function(tp, tn, G0, k, epsilon, beta, N0){
  (P0/(k*tp)+G0*N0*P0/(k+epsilon*P0))/(beta/tn+G0*P0/(k+epsilon*P0))
}
Ss_pre = function(tp, tn, G0, k, epsilon, beta, N0){
  Ns = Ns_pre(tp, tn, G0, k, epsilon, beta, N0)
  (P0/q-Ns/tn)/(G0*(Ns-N0))
}
y_pre = function(tp, tn, G0, k, epsilon, beta, N0){
  Ss = Ss_pre(tp, tn, G0, k, epsilon, beta, N0)
  Ns = Ns_pre(tp, tn, G0, k, epsilon, beta, N0)
  1/tp + 1/tn +G0*Ss/(1+epsilon*Ss) - G0*(Ns-N0)/(1+epsilon*Ss)
}
z_pre = function(tp, tn, G0, k, epsilon, beta, N0){
  Ss = Ss_pre(tp, tn, G0, k, epsilon, beta, N0)
  Ns = Ns_pre(tp, tn, G0, k, epsilon, beta, N0)
  (1/tp - G0*(Ns-N0)/(1+epsilon*Ss))*(1/tn+G0*Ss/(1+epsilon*Ss))+
    (beta/tn+G0*Ss/(1+epsilon*Ss))*G0*(Ns-N0)/(1+epsilon*Ss)
}
# stationary I
I_pre = function(tp, tn, G0, k, epsilon, beta, N0){
  Ss = Ss_pre(tp, tn, G0, k, epsilon, beta, N0)
  Ns = Ns_pre(tp, tn, G0, k, epsilon, beta, N0)
  q*(Ns/tn+G0*(Ns-N0)*P0/(k+epsilon*P0))/ita+Ith
}
# stationary P
P_pre = function(tp, tn, G0, k, epsilon, beta, N0){
  Ns = Ns_pre(tp, tn, G0, k, epsilon, beta, N0)
  k*(P0/q-Ns/tn)/G0/(Ns-N0)
}
Hf_fn = function(tp, tn, G0, k, epsilon, beta, N0, f){
  y = y_pre(tp, tn, G0, k, epsilon, beta, N0)
  z = z_pre(tp, tn, G0, k, epsilon, beta, N0)
  h = z^2/(4*pi^2*f^2*y^2+z^2+(4*pi^2*f^2)^2-8*pi^2*f^2*z)
}
Hf_loss = function(tp, tn, G0, k, epsilon, beta, N0, f, Hf2){
  (Hf_fn(tp, tn, G0, k, epsilon, beta, N0, f)-Hf2)^2
}
loss_fn = function(tp, tn, G0, k, epsilon, beta, N0, f, Hf2){
  I1 = I_pre(tp, tn, G0, k, epsilon, beta, N0)
  P1 = P_pre(tp, tn, G0, k, epsilon, beta, N0)
  Hf_loss(tp, tn, G0, k, epsilon, beta, N0, f, Hf2)+(I1/I-1)^2+(P1/P0-1)^2
}

#### get parameter from python minimize ####
tp=1.29856963e-11
tn=1.00854018e-09
G0=1.63275341e+06
k=1.56080625e-08
epsilon=-2.90219234e-10
beta=-1.96712005e+00
N0=7.96482412e+05

#### theta ####
input2 = read.csv('data/L-I-20C.csv')
abc = as.matrix(data.frame(a = input2$P/1000, b = input2$I/1000, c = input2$I*input2$U/1000-input2$P/1000))
theta = c(3.645557e-01, -3.227650e-03,  4.190533e+03,  0.000000e+00, 
          -2.536997e-05, 2.908151e-07, -8.346502e-10,  9.435902e-13)
pre_fn = function(theta, abc){
  d = T0+(abc[3]+abc[1])*theta[3]
  C0 = theta[1]*(theta[4]+theta[5]*d+theta[6]*d^2+theta[7]*d^3+theta[8]*d^4+theta[2]-abc[2])
  C1 = 1-theta[1]*theta[3]*(theta[5]+2*theta[6]*d+3*theta[7]*d^2+4*theta[8]*d^3)
  C2 = theta[1]*theta[3]^2*(theta[6]+3*theta[7]*d+6*theta[8]*d^2)
  C3 = -theta[1]*theta[3]^3*(theta[7]+4*theta[8]*d)
  C4 = theta[1]*theta[8]*theta[3]^4
  polyroot(c(C0,C1,C2,C3,C4))[1]
}

#### fix I at 7.5mA ####
T0=20+273
I=7.5e-3
ff = f*1e-9
Hf_fn = function(tp, tn, G0, k, epsilon, beta, N0, f){
  y = y_pre(tp, tn, G0, k, epsilon, beta, N0)
  z = z_pre(tp, tn, G0, k, epsilon, beta, N0)
  h = z^2/(4*pi^2*f^2*y^2+z^2+(4*pi^2*f^2)^2-8*pi^2*f^2*z)
  log10(sqrt(h)*h0)*20
}
P0 = pre_fn(theta, abc[abc[,2]==I,])
Hf_pre = Hf_fn(tp, tn, G0, k, epsilon, beta, N0, f)
plot(s21~ff, type='l', xlab = 'f(GHz)', ylab='S21(dB)')
lines(Hf_pre~ff, col =2)
legend('topright', c('Measured', 'Fitted'), col=1:2, lty=1)

# min max
maxhf = max(Re(Hf_pre))
minhf = min(Re(Hf_pre))
for (i in 1:6){
  T0 = i*10+273
  P0 = pre_fn(theta, abc[abc[,2]==I,])
  Hf_pre = Hf_fn(tp, tn, G0, k, epsilon, beta, N0, f)
  maxhf = max(Re(Hf_pre), maxhf)
  minhf = min(Re(Hf_pre), minhf)
}
T0=20+273
P0 = pre_fn(theta, abc[abc[,2]==I,])
Hf_pre = Hf_fn(tp, tn, G0, k, epsilon, beta, N0, f)
plot(s21~ff, type='l', ylim = c(minhf, maxhf), main = 'S21 Curves with I=7.5mA',
     xlab = 'f(GHz)', ylab='S21(dB)')
#lines(Hf_pre~f, col=2)
for (i in 1:6){
  T0 = i*10+273
  P0 = pre_fn(theta, abc[abc[,2]==I,])
  Hf_pre = Hf_fn(tp, tn, G0, k, epsilon, beta, N0, f)
  lines(Hf_pre~ff, col = i+1, lty = i+1)
}
legend('topright', c("Measured",paste(rep('T0 =', 6), (1:6)*10, 'degrees Celsius')), 
       col = 1:7, lty=1:7, cex=1, text.width = 6.6)

#### fix temperature at 20 ####
T0 = 20+273
I=7.5e-3
P0 = abc[abc[,2]==I,1]
Hf_pre = Hf_fn(tp, tn, G0, k, epsilon, beta, N0, f)
maxhf = max(Re(Hf_pre))
minhf = min(Re(Hf_pre))
for (i in 1:9){
  I = 3.5e-3+1e-3*as.double(i)
  P0 = abc[abs(abc[,2]-I)<1e-6,1]
  Hf_pre = Hf_fn(tp, tn, G0, k, epsilon, beta, N0, f)
  maxhf = max(Re(Hf_pre), maxhf)
  minhf = min(Re(Hf_pre), minhf)
}
T0=20+273
I=7.5e-3
P0 = abc[abc[,2]==I,1]
Hf_pre = Hf_fn(tp, tn, G0, k, epsilon, beta, N0, f)
plot(s21~ff, type='l', ylim = c(minhf, maxhf), main = 'S21 Curves with Ta=20 degrees Celsius',
     xlab = 'f(GHz)', ylab='S21(dB)')
#lines(Hf_pre~f, col=2)
for (i in 1:9){
  I = 3.5e-3+1e-3*as.double(i)
  P0 = abc[abs(abc[,2]-I)<1e-6,1]
  Hf_pre = Hf_fn(tp, tn, G0, k, epsilon, beta, N0, f)
  lines(Hf_pre~ff, col = i+1, lty = i+1)
}
legend('topright', c("Measured",paste(rep('I =', 9), (3.5e-3+1e-3*(1:9))*1000, 'mA')), 
       col = 1:10, lty=1:10, cex=1, text.width = 3.5)

