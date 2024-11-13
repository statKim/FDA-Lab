from enum import Enum, auto

import torch
import torch.nn as nn
import torch.nn.functional as F

import numpy as np
import torch.optim as optim



class NoiseType(Enum):
    DIAGONAL = auto() #\Sigma = diag(\sigma^(k))
    ISOTROPIC = auto() #\Sigma = \sigma \times I
    ISOTROPIC_ACROSS_CLUSTERS = auto()
    FIXED = auto()


class MixtureDensityNetwork(nn.Module):
    """
    Mixture density network.

    [ Bishop, 1994 ]

    Parameters
    ----------
    dim_in: int; dimensionality of the covariates
    dim_out: int; dimensionality of the response variable
    n_components: int; number of components in the mixture model
    """
    def __init__(self, dim_in, dim_out, n_components, hidden_dim, noise_type=NoiseType.DIAGONAL, fixed_noise_level=None):
        super().__init__()
        assert (fixed_noise_level is not None) == (noise_type is NoiseType.FIXED)
        num_sigma_channels = {
            NoiseType.DIAGONAL: dim_out * n_components,
            NoiseType.ISOTROPIC: n_components,
            NoiseType.ISOTROPIC_ACROSS_CLUSTERS: 1,
            NoiseType.FIXED: 0,
        }[noise_type]
        self.dim_in, self.dim_out, self.n_components = dim_in, dim_out, n_components
        self.noise_type, self.fixed_noise_level = noise_type, fixed_noise_level
        self.pi_network = nn.Sequential(
            nn.Linear(dim_in, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, n_components),
        )
        self.normal_network = nn.Sequential(
            nn.Linear(dim_in, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, dim_out * n_components + num_sigma_channels)
        )

    def forward(self, x, eps=1e-6):
        #
        # Returns
        # -------
        # log_pi: (bsz, n_components)
        # mu: (bsz, n_components, dim_out)
        # sigma: (bsz, n_components, dim_out)
        #
        log_pi = torch.log_softmax(self.pi_network(x), dim=-1)
        normal_params = self.normal_network(x)
        mu = normal_params[..., :self.dim_out * self.n_components]
        sigma = normal_params[..., self.dim_out * self.n_components:]
        if self.noise_type is NoiseType.DIAGONAL:
            sigma = torch.exp(sigma + eps) # add eps to make it non-zero
        if self.noise_type is NoiseType.ISOTROPIC:
            sigma = torch.exp(sigma + eps).repeat(1, self.dim_out)
        if self.noise_type is NoiseType.ISOTROPIC_ACROSS_CLUSTERS:
            sigma = torch.exp(sigma + eps).repeat(1, self.n_components * self.dim_out)
        if self.noise_type is NoiseType.FIXED:
            sigma = torch.full_like(mu, fill_value=self.fixed_noise_level)
        mu = mu.reshape(-1, self.n_components, self.dim_out)
        sigma = sigma.reshape(-1, self.n_components, self.dim_out)
        return log_pi, mu, sigma

    def loss(self, x, y):
        log_pi, mu, sigma = self.forward(x)
        z_score = (y.unsqueeze(1) - mu) / sigma
        normal_loglik = (
            -0.5 * torch.einsum("bij,bij->bi", z_score, z_score)
            -torch.sum(torch.log(sigma), dim=-1)
        )
        loglik = torch.logsumexp(log_pi + normal_loglik, dim=-1)
        return -loglik

    def sample(self, x):
        log_pi, mu, sigma = self.forward(x)
        cum_pi = torch.cumsum(torch.exp(log_pi), dim=-1)
        rvs = torch.rand(len(x), 1).to(x)
        rand_pi = torch.searchsorted(cum_pi, rvs)
        rand_normal = torch.randn_like(mu) * sigma + mu
        samples = torch.take_along_dim(rand_normal, indices=rand_pi.unsqueeze(-1), dim=1).squeeze(dim=1)
        return samples

    def dens(self, x, y):
        log_pi, mu, sigma = self.forward(x)
        z_score = (y.unsqueeze(1) - mu) / sigma
        normal_loglik = (
            -0.5 * torch.einsum("bij,bij->bi", z_score, z_score)
            -torch.sum(torch.log(sigma), dim=-1)
        )
        loglik = torch.logsumexp(log_pi + normal_loglik, dim=-1)
        density = torch.exp(loglik)
        return density


# Fit Mixture Density Network
def fit_mixture_density_network(x, y, n_components=4, hidden_dim=30, seed=100):
    x = torch.Tensor(x)
    y = torch.Tensor(y)

    torch.manual_seed(seed)
    model = MixtureDensityNetwork(6, 1, n_components=n_components, hidden_dim=hidden_dim, noise_type=NoiseType.DIAGONAL)
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, 2000)

    for i in range(2000):
        optimizer.zero_grad()
        loss = model.loss(x, y).mean()
        loss.backward()
        optimizer.step()
        scheduler.step()
        if i % 200 == 0:
            # logger.info(f"Iter: {i}\t" + f"Loss: {loss.data:.2f}")
            print(loss)
            
    return model


def predict_mixture_density_network(model, x_new, y_new):
    x_new = torch.Tensor(x_new)
    y_new = torch.Tensor(y_new)
    
    f_hat = model.dens(x_new, y_new).detach().numpy()
    
    return f_hat
  
  

