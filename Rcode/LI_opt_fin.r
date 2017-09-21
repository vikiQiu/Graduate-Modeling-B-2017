#### read data ####
input = read.csv('data/L-I-20C.csv')
abc = as.matrix(data.frame(a = input$P/1000, b = input$I/1000, c = input$I*input$U/1000-input$P/1000))
N = nrow(abc)
T0 = 20 + 273

#### calculate x2 ####
# theta0: initial value
theta0 = c(1, -0.000946, 2.6e3+273, 0, -2.545e-5, 2.908e-7,
           -2.531e-10, 1.022e-12)
# start loss
loss0 = apply(abc, 1, loss_fn, theta = theta0)

#### parameter ####
maxIter = 3000
tol = 1e-11
batch_size = N
n_theta = length(theta0)
thetas = matrix(0, maxIter+1, n_theta)
thetas[1,] = theta0
losses = rep(0, maxIter)
losses[1] = mean(loss0)
A_fn = function(theta, abc) {
  d = T0+abc[3]*theta[3]
  theta[1]*(abc[2]-theta[2]-theta[4]-theta[5]*d-theta[6]*d^2
            -theta[7]*d^3-theta[8]*d^4)
}
loss_fn = function(theta, abc){
  A = A_fn(theta, abc)
  (A-abc[1])^2
}
train = function(start, end){
  for (i in start:end) {
    print(paste('Iter',i))
    tmp_theta = rep(0, n_theta)
    for (batch in sample(1:N, batch_size)){
      # print(paste('batch =',batch))
      A_pre = A_fn(thetas[i,], abc[batch,])
      A = abc[batch,1]
      c = abc[batch,3]
      d = T0+c*thetas[i,3]
      tmp_theta[1] = tmp_theta[1] + 2*(A_pre-A)*A_pre/thetas[i,1]
      tmp_theta[2] = tmp_theta[2] - 2*(A_pre-A)*thetas[i,1]
      tmp_theta[3] = tmp_theta[3] - 2*(A_pre-A)*thetas[i,1]*(
        thetas[i,5]*c+2*thetas[i,6]*d*c+3*thetas[i,7]*d^2*c+4*thetas[i,8]*d^3*c)
      for (j in 5:8){
        tmp_theta[j] = tmp_theta[j]-
          2*(A_pre-A)*thetas[i,1]*d^(j-4)
      }
    }
    tmp_theta = tmp_theta/batch_size
    print(tmp_theta)
    thetas[i+1,] = thetas[i,] - learning_rate*tmp_theta
    print(thetas[i+1,])
    losses[i+1] = mean(apply(abc, 1, loss_fn, theta = thetas[i+1,]))
    print(paste('loss =',losses[i+1]))
    if(losses<tol) return(list(thetas=thetas, losses=losses));
  }
  return(list(thetas=thetas, losses=losses))
}

#### learning rate ####
learning_rate = rep(5e-17, n_theta)
learning_rate[8] = 2e-23

#### train ####
res = train(1,50)
ts.plot(res$losses[1:50])
thetas = res$thetas
losses = res$losses
ts.plot(thetas[1:50, 5]) # 1, 3, 4 not update

# 3.28069028538476e-07
learning_rate[1] = 8e4
res = train(50,100)
ts.plot(res$losses[40:100])
thetas = res$thetas
losses = res$losses
ts.plot(thetas[45:100, 7]) # 2, 3, 4 not update

# 1.96792066565275e-08
learning_rate[1] = 3e4
learning_rate[2] = 7e-1
res = train(100,150)
ts.plot(res$losses[80:150])
thetas = res$thetas
losses = res$losses
ts.plot(thetas[95:150, 1])

# 2.43625707009772e-10
learning_rate[1] = 2e4
learning_rate[2] = 5e-1
learning_rate[3] = 3e12
res = train(150,300)
ts.plot(res$losses[145:300])
thetas = res$thetas
losses = res$losses
ts.plot(thetas[145:300, 3])

# 2.14164175440555e-10
learning_rate[1] = 2e4
learning_rate[2] = 5e-1
learning_rate[3] = 6e11
learning_rate[5] = 1e-5
res = train(300,500)
ts.plot(res$losses[250:500])
thetas = res$thetas
losses = res$losses
ts.plot(thetas[200:500, 8])

# 2.1413500835843e-10
learning_rate[1] = 2e4
learning_rate[2] = 5e-1
learning_rate[3] = 6e11
learning_rate[5] = 2e-5
learning_rate[6] = 2e-11
res = train(500,550)
ts.plot(res$losses[490:550])
thetas = res$thetas
losses = res$losses
ts.plot(thetas[450:550, 6])

# 2.13975798009352e-10
learning_rate[1] = 2e4
learning_rate[2] = 5e-1
learning_rate[3] = 6e11
learning_rate[5] = 2e-5
learning_rate[6] = 1e-11
learning_rate[7] = 6e-17
res = train(550,900)
ts.plot(res$losses[500:900])
thetas = res$thetas
losses = res$losses
ts.plot(thetas[450:900, 7])
loss = mean(apply(abc, 1, loss_fn, theta = thetas[900,]))
print(paste('loss=',loss))

# 2.13914103842991e-10
learning_rate[1] = 2e4
learning_rate[2] = 5e-1
learning_rate[3] = 6e11
learning_rate[5] = 2e-5
learning_rate[6] = 1.5e-11
learning_rate[7] = 1e-17
learning_rate[8] = 1e-22
res = train(900,1000)

#### plot ####
par(mfrow = c(2,2))
plot(1:20, res$losses[1:20], type='l', xlab = 'Iteration(1-20)', ylab='Loss', main='(a)')
plot(20:100, res$losses[20:100], type='l', xlab = 'Iteration(20-100)', ylab='Loss', main='(b)')
plot(100:200, res$losses[100:200], type='l', xlab = 'Iteration(100-200)', ylab='Loss', main = '(c)')
plot(200:1000, res$losses[200:1000], type='l', xlab = 'Iteration(200-1000)', ylab='Loss', main = '(d)')
par(mfrow = c(1,1))
thetas = res$thetas
losses = res$losses
ts.plot(thetas[1:1000,2])
loss = mean(apply(abc, 1, loss_fn, theta = thetas[1000,]))
print(paste('loss=',loss))

#### predict ####
theta1 = thetas[1000,]
pre_fn = function(theta, abc){
  d = T0+(abc[3]+abc[1])*theta[3]
  C0 = theta[1]*(theta[4]+theta[5]*d+theta[6]*d^2+theta[7]*d^3+theta[8]*d^4+theta[2]-abc[2])
  C1 = 1-theta[1]*theta[3]*(theta[5]+2*theta[6]*d+3*theta[7]*d^2+4*theta[8]*d^3)
  C2 = theta[1]*theta[3]^2*(theta[6]+3*theta[7]*d+6*theta[8]*d^2)
  C3 = -theta[1]*theta[3]^3*(theta[7]+4*theta[8]*d)
  C4 = theta[1]*theta[8]*theta[3]^4
  polyroot(c(C0,C1,C2,C3,C4))[1]
}
theta_ind = 1000
theta1 = thetas[theta_ind,] #c(3.645557e-01, -3.227650e-03,  4.190533e+03,  0.000000e+00, -2.536997e-05, 2.908151e-07, -8.346502e-10,  9.435902e-13)
T0=10+273
pre2 = apply(abc, 1, pre_fn, theta=theta1)
plot(abc[,2], pre2, type='l', xlab = 'I', ylab = 'P', 
     main = 'L-I Curves Under Different Temperature')
for (i in 2:9){
  T0 = i*10+273
  pre2 = apply(abc, 1, pre_fn, theta=thetas[theta_ind,])
  lines(abc[,2], pre2, col=i, lty=i)
}
legend('topleft', paste(rep('T0 =', 8), (1:9)*10), col = 1:9, lty=1:9)
#write.csv(thetas[1:theta_ind,], 'res/thetas0919_1.csv', row.names = F)
#write.csv(losses[1:theta_ind], 'res/losses0919_1.csv', row.names = F)

