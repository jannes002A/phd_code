function grad_tildeI=get_gradtildeI(Ecost,Ehit,EgradJ,EgradSh)
%This function computes the gradient of the discretized value function,
%'tilde I'

grad_tildeI = EgradJ + (Ecost+Ehit)*EgradSh;

