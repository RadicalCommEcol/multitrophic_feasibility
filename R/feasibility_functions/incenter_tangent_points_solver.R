# Given 3 curves (N*_i=0, N*_j=0, N*_k=0), its incenter I and the radius the incircle (a function
# of the angle theta) we use this functions to obtain Ti, Tj, Tk: the point in N*_i=0, the one in
# N*_j=0 and that of N*_k=0 that are tangent to the (d-1)-spherical cap, centered at I.

incenter_tangent_points_solver <- function(x){ 
  
  solution <- NULL # we store here the resulting values
  
  dimensions <- ncol(A_int)
  
  # Incenter is a global variable: I
  
  # Great circle distancea between I and tangent points = angle theta
  dis_IT1_minus_costheta <- x[1]
  dis_IT2_minus_costheta <- x[2]
  dis_IT1_minus_costheta <- x[3]
  
  # Ti in unit-ball
  T1 <- x[(4+0*dimensions):(3+1*dimensions)]
  T2 <- x[(4+1*dimensions):(3+2*dimensions)]
  T3 <- x[(4+2*dimensions):(3+3*dimensions)]
  
  # # Tang. point in unit-ball
  # solution <- c(solution,sum(T1*T1)-1)
  # solution <- c(solution,sum(T2*T2)-1)
  # solution <- c(solution,sum(T3*T3)-1)
  
  
  # Great circle distancea
  
  I_T1 <- I*T1
  module_T1 <- sqrt(sum(T1*T1))
  I_T1 <- I_T1/module_T1
  
  I_T2 <- I*T2
  module_T2 <- sqrt(sum(T2*T2))
  I_T2 <- I_T2/module_T2
  
  I_T3 <- I*T3
  module_T3 <- sqrt(sum(T3*T3))
  I_T3 <- I_T3/module_T3
  
  solution <- c(solution, sum(I_T1)-cos_theta)
  solution <- c(solution, sum(I_T2)-cos_theta)
  solution <- c(solution, sum(I_T3)-cos_theta)
  
  # Director vectors
  
  dir_vec_1 <- as.numeric(curves_data[1,(3+dimensions):(2+2*dimensions)])
  dir_vec_2 <- as.numeric(curves_data[2,(3+dimensions):(2+2*dimensions)])
  dir_vec_3 <- as.numeric(curves_data[3,(3+dimensions):(2+2*dimensions)])
  
  
  # relations between I and the tang. points
  I_dir_vec_1 <- I - dir_vec_1*sum(I*dir_vec_1)
  module_I_dir_vec_1 <- sqrt(sum(I_dir_vec_1*I_dir_vec_1))
  
  solution <- c(solution,T1[1]-I_dir_vec_1[1]/module_I_dir_vec_1)
  solution <- c(solution,T1[2]-I_dir_vec_1[2]/module_I_dir_vec_1)
  solution <- c(solution,T1[3]-I_dir_vec_1[3]/module_I_dir_vec_1)
  
  I_dir_vec_2 <- I - dir_vec_2*sum(I*dir_vec_2)
  module_I_dir_vec_2 <- sqrt(sum(I_dir_vec_2*I_dir_vec_2))
  
  solution <- c(solution,T2[1]-I_dir_vec_2[1]/module_I_dir_vec_2)
  solution <- c(solution,T2[2]-I_dir_vec_2[2]/module_I_dir_vec_2)
  solution <- c(solution,T2[3]-I_dir_vec_2[3]/module_I_dir_vec_2)
  
  I_dir_vec_3 <- I - dir_vec_3*sum(I*dir_vec_3)
  module_I_dir_vec_3 <- sqrt(sum(I_dir_vec_3*I_dir_vec_3))
  
  solution <- c(solution,T3[1]-I_dir_vec_3[1]/module_I_dir_vec_3)
  solution <- c(solution,T3[2]-I_dir_vec_3[2]/module_I_dir_vec_3)
  solution <- c(solution,T3[3]-I_dir_vec_3[3]/module_I_dir_vec_3)  
  
  return(solution)
  
}
