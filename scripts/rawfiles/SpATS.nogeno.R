##############################################################################################################################################
construct.block <-
function(A1,A2,A3,A4) {
	block <- rbind(cbind(A1,A2), cbind(A3,A4))
	return(block)
}
##############################################################################################################################################
construct.matrices <-
function(X, Z, z, w, GLAM) {
	if(GLAM) {
		XtX. <- XtX(X,w)
		XtZ. <- XtZ(X,Z,w)
		ZtX. <- t(XtZ.)
		ZtZ. <- ZtZ(Z,w)
		Zty. = Zty(Z,z,w)
		Xty. = Xty(X,z,w)
		yty. <- sum((z^2)*w)
		ZtXtZ = rbind(XtZ., ZtZ.)
		u <- c(Xty.,Zty.)
	} else {
		XtW. = t(X*w)
		XtX. = XtW.%*%X
		XtZ. = XtW.%*%Z
		ZtX. = t(XtZ.)
		ZtW. =  t(Z*w)
		ZtZ. = ZtW.%*%Z
		Xty. = XtW.%*%z
		Zty. = ZtW.%*%z
		yty. <- sum((z^2)*w)
		ZtXtZ = rbind(XtZ., ZtZ.)
		u <- c(Xty.,Zty.)
	}
	res <- list(XtX. = XtX., XtZ. = XtZ., ZtX. = ZtX., ZtZ. = ZtZ., Xty. = Xty., Zty. = Zty., yty. = yty., ZtXtZ = ZtXtZ, u = u)
}
##############################################################################################################################################
SpATS.nogeno <- 
function(response, spatial, fixed = NULL, random = NULL, data, family = gaussian(), offset = 0, weights = NULL, control = controlSpATS()) {
	control <- do.call("controlSpATS", control)
	
	if (control$monitoring) start = proc.time()[3]
	
	weights <- as.vector(weights)
	if(is.null(weights)) weights = rep(1, nrow(data))
	if(length(offset) == 1) offset <- rep(offset, nrow(data))

	if(inherits(fixed, "character"))
		fixed <- as.formula(fixed)
	if(inherits(random, "character"))
		random <- as.formula(random)
	if(inherits(spatial, "character"))
		spatial <- as.formula(spatial)
	
	sf <- SpATS:::interpret.SpATS.formula(spatial)
	
	# NAs in the covariates
	model.terms <- c(sf$x.coord, sf$y.coord, if(!is.null(fixed)) attr(terms.formula(fixed), "term.labels"),  if(!is.null(random)) attr(terms.formula(random),"term.labels"))
	na.ind <- apply(is.na(data[,model.terms]), 1, any)
	na.pos <- (1:nrow(data))[!na.ind]
	weights <- weights*(!na.ind)*(!is.na(data[,response]))

	data.na <- data[!na.ind,]
	weights.na <- weights[!na.ind]
	offset.na <-  offset[!na.ind]
	
	y <- data.na[,response]
	nobs <- length(y[weights.na != 0])
	
	###################################################################################################################
	# Construct design matrices
	###################################################################################################################
	MMns <- NULL
    init.var <- gg <- dim <- list()
	int <- rep(1, nrow(data))
    dim.int <- c(Intercept = 1)
    attr(dim.int, "random") <- FALSE
    attr(dim.int, "spatial") <- FALSE
    MMns <- cbind(MMns, Intercept = int)
    dim <- c(dim, list(dim.int)) 
    if (!is.null(fixed)) {
        fixed.part <- SpATS:::construct.fixed.part(formula = fixed, data = data)
        MMns <- cbind(MMns, fixed.part$X)
        dim <- c(dim, list(fixed.part$dim))
    }
    if (!is.null(random)) {
        random.part <- SpATS:::construct.random.part(formula = random, 
            data = data)
        MMns <- cbind(MMns, random.part$Z)
        dim <- c(dim, list(random.part$dim))
        gg <- c(gg, list(random.part$g))
        init.var <- c(init.var, list(rep(0.001, length(random.part$init.var))))
    }
    spat.part <- SpATS:::construct.2d.pspline(formula = spatial, data = data)
    MMns <- cbind(MMns, spat.part$X, spat.part$Z)
    dim <- c(dim, list(spat.part$dim$fixed), list(spat.part$dim$random))
    gg <- c(gg, list(spat.part$g))
    init.var <- c(init.var, list(spat.part$init.var))
    
    g <- SpATS:::construct.capital.lambda(gg)
    random. <- unlist(lapply(dim, SpATS:::get.attribute, "random"))
    spatial. <- unlist(lapply(dim, SpATS:::get.attribute, "spatial"))
    spat.part$terms$fixed$pos <- SpATS:::create.position.indicator(unlist(dim),!random. & spatial.)
    spat.part$terms$random$pos <- SpATS:::create.position.indicator(unlist(dim),random. & spatial.)
    res <- list()
    res$fixed.part <- if (!is.null(fixed)){
        fixed.part
    } else{
    	NULL
    }
    res$random.part <- if (!is.null(random)){
        random.part
    } else {
    	NULL
    }
    res$spat.part <- spat.part
    res$terms$spatial <- spat.part$terms
    res$terms$fixed <- if (!is.null(fixed)){ 
        fixed.part$terms
    } else {
    	NULL
    }
    res$terms$random <- if (!is.null(random)){ 
        random.part$terms
    } else {
    	NULL
    }
    res$g <- g
    res$dim <- dim
    res$init.var <- init.var
    res$MM <- list(MMns = MMns)
	###################################################################################################################
	###################################################################################################################
	MM <- res
	ldim <- MM$dim
	random. <- unlist(lapply(ldim, SpATS:::get.attribute, "random"))
	spatial. <- unlist(lapply(ldim, SpATS:::get.attribute, "spatial"))
	dim <- unlist(ldim)
	g <- MM$g

	# Nominal dimension
	dim.nom <- SpATS:::obtain.nominal.dimension(MM$MM$MMns, dim, random., spatial., weights.na)
	
	random.pos <- SpATS:::create.position.indicator(dim, random.)
	fixed.pos <- SpATS:::create.position.indicator(dim, !random.) 	
	df.fixed <- sum(dim[!random.])
	
	X <- MMns[,fixed.pos,drop = FALSE]
	Z <- MMns[,random.pos,drop = FALSE]
	
	np <- c(ncol(X), ncol(Z))
	D <- diag(c(rep(0,np[1]), rep(1,sum(np[2]))))
	# Fixed and random effects
	# bold = rep(0, sum(dim, na.rm = TRUE))
	# Fit the model
	# eta <- MMns%*%bold + offset.na
	# mu <- family$linkinv(eta)
	mustart <- etastart <- NULL
	eval(family$initialize)
	mu <- mustart
	eta <- family$linkfun(mustart)
	
	# Initialize variance components
	la <- c(1, unlist(MM$init.var))
	# Initialize deviance and psi
	devold <- 1e10
	psi <- la[1]
	
	if(control$monitoring > 1) {
		cat("Effective dimensions\n")
		cat("-------------------------\n")
		cat(sprintf("%1$3s %2$12s","It.","Deviance"), sep = "")
		cat(sprintf("%12s", names(g)), sep = "")
		cat("\n")
	}
	for (iter in 1:control$maxit) {
		deriv <- family$mu.eta(eta)
		z <- (eta - offset.na) + (y - mu)/deriv
		w <- as.vector(deriv^2/family$variance(mu))
		w <- w*weights.na
		z[!weights.na] <- 0
		mat <- construct.matrices(X, Z, z, w, FALSE)
		
		if(control$monitoring) start1 <- proc.time()[3]
		for (it in 1:control$maxit) {
			# Build penalty matrix: block diagonal matrix
			Ginv <- vector(length = sum(dim[random.], na.rm = TRUE))
			for (i in 1:length(g)) {
				Ginv <- Ginv + (1/la[i+1])*g[[i]]
			}
			G = 1/Ginv
			V <- construct.block(mat$XtX., t(mat$ZtX.*G), mat$ZtX., t(mat$ZtZ.*G))
		
			H <- (1/la[1])*V + D
			Hinv <- solve(H)
			b.aux <- (1/la[1])*Hinv%*%mat$u
			
			b.fixed <- b.aux[1:np[1]]
			b.random <- G*b.aux[-(1:np[1])]
			
			b <- rep(0, sum(np))
			b[fixed.pos] <- b.fixed
			b[random.pos] <- b.random
			
			dZtPZ <- 1/la[1]*apply((t(Hinv[-(1:np[1]),])*mat$ZtXtZ),2,sum)
			
			ed <- tau <- vector(mode="list")
			for (i in 1:length(g)) {
				g.inv.d <- (1/la[i+1])*g[[i]]
				ed[[i]] <- sum(dZtPZ*(g.inv.d*G^2))
				ed[[i]] <- ifelse(ed[[i]] == 0, 1e-50, ed[[i]])
				tau[[i]] <- sum(b.random^2*g[[i]])/ed[[i]]
				tau[[i]] <- ifelse(tau[[i]] == 0, 1e-50, tau[[i]])
			}
			ssr = mat$yty. - t(c(b.fixed, b.random))%*%(2*mat$u - V%*%b.aux)
			# Compute deviance
			dev <- determinant(H)$modulus + sum(log(la[1]*1/w[w != 0])) + ssr/la[1] + sum(b.random^2*Ginv)
			psinew <- as.numeric((ssr/(nobs - sum(unlist(ed)) - df.fixed)))
			if(family$family == "gaussian" | control$update.psi) {
				psi2 <- psinew
			} else {
				psi2 <- 1
			}
			# New variance components and convergence check
			lanew <- c(psi2, unlist(tau))
			dla = abs(devold - dev)
			if(control$monitoring > 1) {
				cat(sprintf("%1$3d %2$12.6f", it, dev), sep = "")
				cat(sprintf("%12.3f", unlist(ed)), sep = "")
				cat('\n')
			}
			if (dla < control$tolerance) break
			la = lanew
			psi = psinew
			devold = dev
		}
		if (control$monitoring) {
			end1 <- proc.time()[3]
			cat("Timings:\nSpATS", (end1-start1), "seconds\n")
		}
		eta.old <- eta
		eta <- MMns%*%b + offset.na
		mu <- family$linkinv(eta)
		# Convergence criterion: linear predictor
		tol <- sum((eta - eta.old)^2)/sum(eta^2)
		if (tol < control$tolerance | (family$family == "gaussian" & family$link== "identity")) break
	}	
	var.comp <- la[-1]
	eff.dim <- unlist(ed)
	names(var.comp) <- names(eff.dim) <- names(g)
	
	# Effective dimension (fixed + random)
	eff.dim <- c(dim[!random.], eff.dim)
	
	attr(dim, "random") <- random.
	attr(dim, "spatial") <- spatial.
	
	fitted <- rep(NA, nrow(data))
	fitted[!na.ind] <- mu
	
	# Deviance residuals
	dev.residuals <- family$dev.resids(data[,response], fitted, weights)
    s <- attr(dev.residuals, "sign")
    if (is.null(s)) 
        s <- sign(data[,response] - fitted)
    dev.residuals <- sqrt(pmax(dev.residuals, 0))*s
    
	if (control$monitoring) {
		end <- proc.time()[3]
		cat("All process", (end - start), "seconds\n")
	}
	res <- list()
	res$call <- match.call()
	res$data <- cbind(data, weights = weights)
	res$model <- list(response = response, spatial = spatial, fixed = fixed, random = random)
	res$fitted <- fitted
	res$residuals <- dev.residuals
	res$psi <- c(la[1], psi)
	res$var.comp <- var.comp
	res$eff.dim <- eff.dim
	res$dim <- dim
	res$dim.nom <- dim.nom
	res$nobs <- nobs
	res$deviance <- dev
	res$coeff <- b
	random.coeff <- rep(FALSE, length(b))
	random.coeff[SpATS:::create.position.indicator(dim, random.)] <- TRUE
	attr(res$coeff, "random") <- random.coeff	
	# Terms
	res$terms <- MM$terms
	class(res) <- "SpATS"
	res
}
