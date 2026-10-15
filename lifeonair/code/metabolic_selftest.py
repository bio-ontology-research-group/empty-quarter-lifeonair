#!/usr/bin/env python3
"""Small synthetic LPs verify exchange direction and stoichiometric scaling."""
import cobra
from metabolic_simulate import set_medium,solve
for coef in [-1,1,-2,2]:
 m=cobra.Model('exchange_direction');m.solver='glpk'
 a=cobra.Metabolite('cpd00027_e0',compartment='e0')
 b=cobra.Metabolite('cpd00027_c0',compartment='c0')
 ex=cobra.Reaction('EX');ex.add_metabolites({a:coef})
 tr=cobra.Reaction('transport');tr.add_metabolites({a:-1,b:1})
 demand=cobra.Reaction('test_sink');demand.add_metabolites({b:-1})
 m.add_reactions([ex,tr,demand]);m.objective=demand
 for bound in [0,.01,2]:
  set_medium(m,[ex],{'cpd00027':bound});status,val=solve(m)
  assert status=='optimal' and abs(val-bound)<1e-8,(coef,bound,status,val)
 print('PASS exchange coefficient',coef,'zero and bounded uptake')
print('12 exchange sign/scaling tests passed')
