################################################################################################
# Fitting a penalized Bezier curve fitting with general order
# Maintainer: Jae-Kyung Shin
# email: tlsworud60@korea.ac.kr
################################################################################################

# fitting a penalized Bezier curve fitting to the data
# This function proceeds with tuning parameter sequence
penalized_bezier = function(tvec, y, control_points = NULL, dimension = 3,
                            jump_eps = 1e-02, lambdas = 0, step_size = 1, maxiter = 1000,
                            epsilon_iter = 1e-06, verbose = FALSE)
{
   sample_size = length(tvec)
   if (is.null(control_points))
      control_points = y[round(seq(1, sample_size, length = dimension)), ]
   dimension = nrow(control_points)
   degree = dimension - 1
   number_penalty = dimension - 2
   distance_R = rep(0, sample_size)
   # initialization
   beta_t = general_bezier(tvec, control_points)
   for (i in 1 : sample_size)
      distance_R[i] = Acos(dot(beta_t[i, ], y[i, ]))^2
   Rlambda = mean(distance_R)
   if (number_penalty > 0)
      Rlambda = Rlambda + lambdas[1] * sum(rowSums(jump_linear(control_points, 1 : dimension)^2))
   R_step = Rlambda
   Rlambda_stored = Inf
   temp_cp = control_points
   # for information of each lambda
   fit = list()
   number_lambdas = length(lambdas)
   R_list = rep(0, number_lambdas)
   bic_list = rep(0, number_lambdas)
   dimension_list = rep(0, number_lambdas)
   # cda loop
   for (lambda_index in 1 : number_lambdas)
   {
      if (verbose)
         cat(lambda_index, "th lambdas runs \n")
      lambda = lambdas[lambda_index]
      for (iter in 1 : maxiter)
      {
         if (verbose)
            cat(iter, "th iteration runs \n")
         for (j in 1 : dimension)
         {
            gradient = Rgradient_loss(y, tvec, control_points, j)
            if (number_penalty > 0 & lambda > 0)
               gradient = gradient + lambda * R_gradient_penalty(control_points, 1 : dimension, j)$R_grad
            R_step = Rlambda
            step = step_size
            for (iter_step in 1 : 100)
            {
               control_point_tmp = Exp(control_points[j, ], -gradient * step)
               # solve numerical issues
               temp_cp[j, ] = control_point_tmp / norm2(control_point_tmp)
               beta_t = general_bezier(tvec, temp_cp)
               for (i in 1 : sample_size)
                  distance_R[i] = Acos(dot(beta_t[i, ], y[i, ]))^2
               Rlambda = mean(distance_R)
               if (number_penalty > 0 & lambda > 0)
                  Rlambda = Rlambda + lambda * sum(rowSums(jump_linear(temp_cp, 1 : dimension)^2))
               if (R_step >= Rlambda)
                  break
               step = step / 2
            }
            control_points[j, ] = temp_cp[j, ]
         }
         # pruning step
         if (number_penalty > 0 & lambda > 0)
         {
            jump_size = rep(0, number_penalty)
            jump_vector = jump_linear(control_points, 1 : dimension)
            for (k in 1 : number_penalty)
               jump_size[k] = norm2(jump_vector[k, ])^2
            penalty_check = jump_size < jump_eps
            if (sum(penalty_check) > 0)
            {
               prune_index = which(penalty_check)
               control_points = control_points[-(prune_index + 1), ]
               dimension = nrow(control_points)
               number_penalty = dimension - 2
               temp_cp = control_points
               control_points = y[round(seq(1, sample_size, length = dimension)), ]
            }
         }
         # calcuate fitted values
         beta_t = general_bezier(tvec, control_points)
         # copute the squared error
         for (i in 1 : sample_size)
            distance_R[i] = Acos(dot(beta_t[i, ], y[i, ]))^2
         Rlambda = mean(distance_R)
         if (number_penalty > 0 & lambda > 0)
            Rlambda = Rlambda + lambda * sum(rowSums(jump_linear(control_points, 1 : dimension)^2))
         if (verbose)
            cat("R = ", Rlambda, "\n")
         if (abs(Rlambda - Rlambda_stored) < epsilon_iter)
            break
         Rlambda_stored = Rlambda
      }
      for (i in 1 : sample_size)
         distance_R[i] = Acos(dot(beta_t[i, ], y[i, ]))^2
      R = mean(distance_R)
      fit[[lambda_index]] = list(beta_t = beta_t, step_size = step_size,
                                 control_points = control_points)
      bic_list[lambda_index] = sample_size * log(R) + 3 * dimension * log(sample_size)
      # test
      for (i in 1 : sample_size)
         distance_R[i] = Acos(dot(beta_t[i, ], y[i, ]))
      R = mean(distance_R)
      dimension_list[lambda_index] = dimension
      R_list[lambda_index] = R
   }
   fit$bic_list = bic_list
   fit$dimension_list = dimension_list
   fit$R_list = R_list
   return(fit)
}

# This makes the shperical Bezier curves in a recursive way
general_bezier = function(t, control_points)
{
   # initial set-up
   degree = nrow(control_points) - 1
   beziers = list()
   for (j in 1 : degree)
      beziers[[j]] = Geodesic(t, control_points[j, ], control_points[j + 1, ], 0, 1)
   while (degree > 1)
   {
      for (j in 1 : (degree - 1))
         beziers[[j]] = order_bezier(t, beziers[[j]], beziers[[j + 1]])
      degree = degree - 1
   }
   return(beziers[[1]])
}

# (sin(c(t) * (1 - t)) * p_t + sin(c(t) * t) q_t) / sin(c(t))
order_bezier = function(t, p_t, q_t)
{
   c_t = calculate_c_t(p_t, q_t)
   R_1_t_c = calculate_R_s_c_t(c_t, 1 - t)
   R_t_c = calculate_R_s_c_t(c_t, t)
   return(R_1_t_c * p_t + R_t_c * q_t)
}

# This makes the control polygon for a given set of control points
polygon_bezier = function(t, control_points)
{
   gamma = matrix(nrow = 0, ncol = 3)
   for (j in 1 : (nrow(control_points) - 1))
   {
      piece_gamma = Geodesic(t, control_points[j, ], control_points[j + 1, ],
                             0, 1)
      gamma = rbind(gamma, piece_gamma)
   }
   return(gamma)
}

# R(theta, s) = sin(s theta) / sin(theta)
calculate_R_s = function(theta, s)
{
   if (theta == 0)
      return(0)
   else
      return(sin(theta * s) / sin(theta))
}

# calculate Q_s(theta)
calculate_Q_s = function(theta, s)
{
   value_Qs = (sin(s * theta) * cos(theta) - s * cos(s * theta) * sin(theta)) / sin(theta)^3
   return(value_Qs)
}

# calculate c(t) = dist(p(t), q(t))
calculate_c_t = function(p_t, q_t)
{
   c_t = rep(0, nrow(p_t))
   for (i in 1 : nrow(p_t))
      c_t[i] = dist(p_t[i, ], q_t[i, ])
   return(c_t)
}

# calculate R_s(c(t))
calculate_R_s_c_t = function(c_t, s)
{
   R_s_c_t = rep(0, length(s))
   for (i in 1 : length(s))
      R_s_c_t[i] = calculate_R_s(c_t[i], s[i])
   return(R_s_c_t)
}

# psi / sin(psi)
calculate_Apsi = function(psi)
{
   if (sum(psi^2) == 0)
      return(0)
   else
      return(psi / sin(psi))
}

# projection y on sphere onto tangent plane at p
calculate_projection_p = function(p, y)
{
   proj_y = (y - p * dot(p, y)) #/ norm2(y - p * dot(p, y))
   return(proj_y)
}

# This computes the gradient of the linear Bezier curve at a point
gradient_linear_point = function(t, control_points, which, index)
{
   if (sum(which == index) == 0)
      return(matrix(0, 3, 3))
   grad_linear = matrix(0, 3, 3)
   theta = dist(control_points[1, ], control_points[2, ])
   Q_t = calculate_Q_s(theta, t)
   Q_1_t = calculate_Q_s(theta, 1 - t)
   if (index == which[1])
   {
      R_1_t = calculate_R_s(theta, 1 - t)
      grad_linear = R_1_t * diag(1, 3) + Q_1_t * outer(control_points[2, ], control_points[1, ]) +
         Q_t * outer(control_points[2, ], control_points[2, ]) ## Q위치
   }
   else if (index == which[2])
   {
      R_t = calculate_R_s(theta, t)
      grad_linear = R_t * diag(1, 3) + Q_1_t * outer(control_points[1, ], control_points[1, ]) +
         Q_t * outer(control_points[1, ], control_points[2, ]) ## Q위치
   }
   return(grad_linear)
}

# This computes the gradient of the general Bezier curve at a point
gradient_bezier_point = function(t, control_points, index)
{
   # browser()
   degree = nrow(control_points) - 1
   gradients_bezier = list()
   for (j in 1 : degree)
      gradients_bezier[[j]] = gradient_linear_point(t, control_points[j : (j + 1), ], j : (j + 1), index)
   i = 1
   while (degree > 1)
   {
      # browser()
      temp_gradients_bezier = gradients_bezier
      for (j in 1 : degree)
         temp_gradients_bezier[[j]] = t(gradients_bezier[[j]])
      # make bezier curve at a point t
      bezier = list()
      for (j in 1 : degree)
      {
         whichs = j : (j + i)
         bezier[[j]] = general_bezier(t, control_points[whichs, ])
      }
      i = i + 1
      # calculate the gradients of current bezier
      for (j in 1 : (degree - 1))
      {
         theta = dist(bezier[[j]], bezier[[j + 1]])
         R_t = calculate_R_s(theta, t)
         R_1_t = calculate_R_s(theta, 1 - t)
         Q_t = calculate_Q_s(theta, t)
         Q_1_t = calculate_Q_s(theta, 1 - t)
         # browser()
         grad_R_t = Q_t * (t(temp_gradients_bezier[[j]]) %*% t(bezier[[j + 1]]) +
                           t(temp_gradients_bezier[[j + 1]]) %*% t(bezier[[j]]))
         grad_R_1_t = Q_1_t * (t(temp_gradients_bezier[[j]]) %*% t(bezier[[j + 1]]) +
                               t(temp_gradients_bezier[[j + 1]]) %*% t(bezier[[j]]))
         gradients_bezier[[j]] = R_1_t * t(temp_gradients_bezier[[j]]) + R_t * t(temp_gradients_bezier[[j + 1]]) +
                                 grad_R_1_t %*% bezier[[j]] + grad_R_t %*% bezier[[j + 1]]
      }
      degree = degree - 1
   }
   return(t(gradients_bezier[[1]]))
}

# This computes the Riemannian gradient of the loss function at a point
Rgradient_loss_point = function(y, t, control_points, index)
{
   # browser()
   grad_bezier = gradient_bezier_point(t, control_points, index)
   bezier_t = general_bezier(t, control_points)
   phi = dist(y, bezier_t)
   proj = calculate_projection_p(control_points[index, ], t(grad_bezier) %*% y)
   Aphi = calculate_Apsi(phi)
   Rgradient_bezier = - Aphi * proj
   return(Rgradient_bezier)
}

# This computes the Riemannian gradient of the loss function
Rgradient_loss = function(y, t, control_points, index)
{
   # browser()
   N = length(t)
   grad_bezier = matrix(0, 1, 3)
   for (n in 1 : N)
      grad_bezier = grad_bezier + t(Rgradient_loss_point(y[n, ], t[n], control_points, index))
   return(grad_bezier)
}

# This computes the Riemannian gradient of the penalty function
R_gradient_penalty = function(control_points, knots, index)
{
   # initial setting
   dimension = nrow(control_points)
   theta = rep(0, 2)
   delta = rep(0, 2)
   R_gradients = matrix(0, nrow = 3, ncol = 1)
   if (index < dimension - 1)
   {
      for (k in 1 : 2)
      {
         delta[k] = knots[index + k] - knots[index + k - 1]
         theta[k] = dist(control_points[index + k, ], control_points[index + k - 1, ])
      }
      # calculate the constant terms of penalty function
      a = theta[1] / (sin(theta[1]) * delta[1])
      b1 = cos(theta[1]) * theta[1] / (sin(theta[1]) * delta[1])
      b2 = cos(theta[2]) * theta[2] / (sin(theta[2]) * delta[2])
      c = theta[2] / (sin(theta[2]) * delta[2])
      # difference of the corresponding penalty
      d = a * control_points[index, ] - (b1 + b2) * control_points[index + 1, ] + c * control_points[index + 2, ]
      # values of derivetives
      prime1 = (theta[1] * cos(theta[1]) - sin(theta[1])) / (sin(theta[1]))^3 * outer(control_points[index + 1, ], control_points[index, ])
      prime2 = theta[1] / sin(theta[1]) * diag(3)
      # calculate the gradients
      grad = (prime1 + prime2) %*% d
      grad = grad / (norm2(d) * delta[1])
      R_gradients = R_gradients + grad
   }
   if (index > 1 & index < dimension)
   {
      for (k in 1 : 2)
      {
         delta[k] = knots[index + k - 1] - knots[index + k - 2]
         theta[k] = dist(control_points[index + k - 1, ], control_points[index + k - 2 , ])
      }
      # calculate the constant terms of penalty function
      a = theta[1] / (sin(theta[1]) * delta[1])
      b1 = cos(theta[1]) * theta[1] / (sin(theta[1]) * delta[1])
      b2 = cos(theta[2]) * theta[2] / (sin(theta[2]) * delta[2])
      c = theta[2] / (sin(theta[2]) * delta[2])
      # difference of the corresponding penalty
      d = a * control_points[index - 1, ] - (b1 + b2) * control_points[index, ] + c * control_points[index + 1, ]
      # values of derivetives
      prime1 = (theta[2] * cos(theta[2]) - sin(theta[2])) / (sin(theta[2]))^3 * outer(control_points[index + 1, ], control_points[index + 1, ]) / delta[2]
      prime2 = (theta[1] * cos(theta[1]) - sin(theta[1])) / (sin(theta[1]))^3 * outer(control_points[index - 1, ], control_points[index - 1, ]) / delta[1]
      prime5 = theta[2] * cos(theta[2]) / sin(theta[2]) * diag(3) / delta[2]
      prime6 = theta[1] * cos(theta[1]) / sin(theta[1]) * diag(3) / delta[1]
      # calculate the gradients
      grad = (prime1 + prime2 - prime5 - prime6) %*% d
      grad = grad / norm2(d)
      R_gradients = R_gradients + grad
   }
   if (index > 2)
   {
      for (k in 1 : 2)
      {
         delta[k] = knots[index + k - 2] - knots[index + k - 3]
         theta[k] = dist(control_points[index + k - 2, ], control_points[index + k - 3, ])
      }
      # calculate the constant terms of penalty function
      a = theta[1] / (sin(theta[1]) * delta[1])
      b1 = cos(theta[1]) * theta[1] / (sin(theta[1]) * delta[1])
      b2 = cos(theta[2]) * theta[2] / (sin(theta[2]) * delta[2])
      c = theta[2] / (sin(theta[2]) * delta[2])
      # difference of the corresponding penalty
      d = a * control_points[index - 2, ] - (b1 + b2) * control_points[index - 1, ] + c * control_points[index, ]
      # values of derivetives
      prime1 = (theta[2] * cos(theta[2]) - sin(theta[2])) / (sin(theta[2]))^3 * outer(control_points[index - 1, ], control_points[index, ])
      prime2 = theta[2] / sin(theta[2]) * diag(3)
      # calculate the gradients
      grad = (prime1 + prime2) %*% d
      grad = grad / (norm2(d) * delta[2])
      R_gradients = R_gradients + grad
   }
   # calculate the Rimannian gradient and Exp(Rimannian gradient)
   R_gradients = t(as.matrix(calculate_projection_p(control_points[index, ], R_gradients)))
   Exp_R_gradients = Exp(control_points[index, ], - R_gradients) / norm2(Exp(control_points[index, ], - R_gradients))
   return(list(R_grad = R_gradients, Exp_Rg = Exp_R_gradients))
}

# This computes value of jump size for the penalty
jump_linear = function(control_points, knots)
{
   number_penalty = nrow(control_points) - 2
   jump_vector = matrix(0, number_penalty, 3)
   for (j in 1 : number_penalty)
   {
      theta1 = Acos(dot(control_points[j, ], control_points[j + 1, ]))
      theta2 = Acos(dot(control_points[j + 1, ], control_points[j + 2, ]))
      a = theta2 / sin(theta2) / (knots[j + 2] - knots[j + 1])
      c = theta1 / sin(theta1) / (knots[j + 1] - knots[j])
      b1 = cos(theta2) * a
      b2 = cos(theta1) * c
      jump_vector[j, ] = a * control_points[j + 2, ] -
         (b1 + b2) * control_points[j + 1, ] + c * control_points[j, ]
   }
   return(jump_vector)
}
