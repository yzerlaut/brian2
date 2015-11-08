# TODO: handle correctly things like ((a*2)*3)*4 -> a*24
# we can do this by checking as part of Mult handling that all subnodes are mults or divs, ...

import ast
from collections import OrderedDict
import copy
import itertools

from brian2.core.functions import DEFAULT_FUNCTIONS, DEFAULT_CONSTANTS
from brian2.core.variables import Variable
from brian2.parsing.bast import (brian_ast, BrianASTRenderer, dtype_hierarchy, is_boolean_dtype,
                                 brian_dtype_from_dtype)
from brian2.parsing.rendering import NodeRenderer
from brian2.utils.stringtools import get_identifiers, word_substitute

from .statements import Statement

defaults_ns = dict((k, v.pyfunc) for k, v in DEFAULT_FUNCTIONS.iteritems())
defaults_ns.update(**dict((k, v.value) for k, v in DEFAULT_CONSTANTS.iteritems()))


__all__ = ['optimise_statements', 'ArithmeticSimplifier', 'Simplifier']


def evaluate_expr(expr, ns=None):
    if ns is None:
        ns = defaults_ns
    try:
        val = eval(expr, ns)
        return val, True
    except NameError:
        return expr, False


def optimise_statements(scalar_statements, vector_statements, variables):
    boolvars = dict((k, v) for k, v in variables.iteritems()
                    if hasattr(v, 'dtype') and brian_dtype_from_dtype(v.dtype)=='boolean')
    simplifier = Simplifier(variables, scalar_statements)
    new_vector_statements = []
    for stmt in vector_statements:
        new_expr = simplifier.render_expr(stmt.expr)
        new_stmt = Statement(stmt.var, stmt.op, new_expr, stmt.comment,
                             dtype=stmt.dtype,
                             constant=stmt.constant,
                             subexpression=stmt.subexpression,
                             scalar=stmt.scalar)
        #complexity_std = expression_complexity(expr_std)
        idents = get_identifiers(new_expr)
        used_boolvars = [var for var in boolvars.iterkeys() if var in idents]
        if len(used_boolvars):
            bool_space = [[False, True] for var in used_boolvars]
            expanded_expressions = {}
            #complexities = {}
            for bool_vals in itertools.product(*bool_space):
                subs = dict((var, str(val)) for var, val in zip(used_boolvars, bool_vals))
                curexpr = word_substitute(new_expr, subs)
                curexpr = simplifier.render_expr(curexpr)
                key = tuple((var, val) for var, val in zip(used_boolvars, bool_vals))
                expanded_expressions[key] = curexpr
                #complexities[key] = expression_complexity(curexpr)
                # print ', '.join('%s=%s'%(k, v) for k, v in key)
                # print '-> ', curexpr
            new_stmt.used_boolean_variables = used_boolvars
            new_stmt.boolean_simplified_expressions = expanded_expressions
        new_vector_statements.append(new_stmt)
    new_scalar_statements = copy.copy(scalar_statements)
    for expr, name in simplifier.loop_invariants.iteritems():
        dtype_name = simplifier.loop_invariant_dtypes[name]
        if dtype_name=='boolean':
            dtype = bool
        elif dtype_name=='integer':
            dtype = int
        else:
            dtype = float
        new_stmt = Statement(name, ':=', expr, '',
                             dtype=dtype,
                             constant=True,
                             subexpression=False,
                             scalar=True)
        new_scalar_statements.append(new_stmt)
    return new_scalar_statements, new_vector_statements


def create_assumptions_namespace(assumptions):
    ns = defaults_ns.copy()
    for assumption in assumptions:
        try:
            exec assumption in ns
        except NameError:
            pass
    return ns


class ArithmeticSimplifier(BrianASTRenderer):
    def __init__(self, variables, assumptions=None):
        BrianASTRenderer.__init__(self, variables)
        if assumptions is None:
            assumptions = []
        self.assumptions = assumptions
        self.assumptions_ns = create_assumptions_namespace(assumptions)

    def render_node(self, node):
        '''
        Assumes that the node has already been fully processed by BrianASTRenderer
        '''
        node = super(ArithmeticSimplifier, self).render_node(node)
        # can't evaluate vector expressions, so abandon in this case
        if not node.scalar:
            return node
        # try fully evaluating using assumptions
        expr = NodeRenderer().render_node(node)
        val, evaluated = evaluate_expr(expr, self.assumptions_ns)
        if evaluated:
            if node.dtype=='boolean':
                val = bool(val)
                if hasattr(ast, 'NameConstant'):
                    newnode = ast.NameConstant(val)
                else:
                    # None is the expression context, we don't use it so we just set to None
                    newnode = ast.Name(repr(val), None)
            elif node.dtype=='integer':
                val = int(val)
            else:
                val = float(val)
            if node.dtype!='boolean':
                newnode = ast.Num(val)
            newnode.dtype = node.dtype
            newnode.scalar = True
            newnode.complexity = 0
            return newnode
        return node

    def render_BinOp(self, node):
        node.left = self.render_node(node.left)
        node.right = self.render_node(node.right)
        node = super(ArithmeticSimplifier, self).render_BinOp(node)
        left = node.left
        right = node.right
        op = node.op
        # Handle multiplication by 0 or 1
        if op.__class__.__name__=='Mult':
            if left.__class__.__name__=='Num':
                if left.n==0:
                    # must not change the dtype of the output, e.g. handle 0*float->0.0 and 0.0*int->0.0
                    left.dtype = node.dtype
                    if node.dtype=='integer':
                        left.n = 0
                    else:
                        left.n = 0.0
                    return left
                if left.n==1:
                    # only simplify this if the type wouldn't be cast by the operation
                    if dtype_hierarchy[left.dtype]<=dtype_hierarchy[right.dtype]:
                        return right
            if right.__class__.__name__=='Num':
                if right.n==0:
                    # must not change the dtype of the output, e.g. handle 0*float->0.0 and 0.0*int->0.0
                    right.dtype = right.dtype
                    if node.dtype=='integer':
                        right.n = 0
                    else:
                        right.n = 0.0
                    return right
                if right.n==1:
                    # only simplify this if the type wouldn't be cast by the operation
                    if dtype_hierarchy[right.dtype]<=dtype_hierarchy[left.dtype]:
                        return left
        # Handle division by 1, or 0/x
        if op.__class__.__name__=='Div':
            if left.__class__.__name__=='Num':
                if left.n==0:
                    # must not change the dtype of the output, e.g. handle 0/float->0.0 and 0.0/int->0.0
                    left.dtype = node.dtype
                    if node.dtype=='integer':
                        left.n = 0
                    else:
                        left.n = 0.0
                    return left
            if right.__class__.__name__=='Num':
                if right.n==1:
                    # only simplify this if the type wouldn't be cast by the operation
                    if dtype_hierarchy[right.dtype]<=dtype_hierarchy[left.dtype]:
                        return left
        # Handle addition of 0
        if op.__class__.__name__=='Add':
            if left.__class__.__name__=='Num':
                if left.n==0:
                    # only simplify this if the type wouldn't be cast by the operation
                    if dtype_hierarchy[left.dtype]<=dtype_hierarchy[right.dtype]:
                        return right
            if right.__class__.__name__=='Num':
                if right.n==0:
                    # only simplify this if the type wouldn't be cast by the operation
                    if dtype_hierarchy[right.dtype]<=dtype_hierarchy[left.dtype]:
                        return left
        # Handle subtraction of 0
        if op.__class__.__name__=='Sub':
            if right.__class__.__name__=='Num':
                if right.n==0:
                    # only simplify this if the type wouldn't be cast by the operation
                    if dtype_hierarchy[right.dtype]<=dtype_hierarchy[left.dtype]:
                        return left
        # simplify e.g. 2*float to 2.0*float to make things more explicit: not strictly necessary
        # but might be useful for some codegen targets
        if node.dtype=='float' and op.__class__.__name__ in ['Mult', 'Add', 'Sub', 'Div']:
            for subnode in [node.left, node.right]:
                if subnode.__class__.__name__=='Num':
                    subnode.dtype = 'float'
                    subnode.n = float(subnode.n)
        return node


class Simplifier(BrianASTRenderer):
    def __init__(self, variables, scalar_statements):
        BrianASTRenderer.__init__(self, variables)
        self.loop_invariants = OrderedDict()
        self.loop_invariant_dtypes = {}
        self.n = 0
        self.node_renderer = NodeRenderer(use_vectorisation_idx=False)
        self.arithmetic_simplifier = ArithmeticSimplifier(variables)
        self.scalar_statements = scalar_statements

    def render_expr_with_additional_assumptions(self, expr, additional_assumptions):
        cur_arithmetic_simplifier = self.arithmetic_simplifier
        self.arithmetic_simplifier = ArithmeticSimplifier(self.variables,
                                                          assumptions=additional_assumptions)
        expr = self.render_expr(expr)
        self.arithmetic_simplifier = cur_arithmetic_simplifier
        return expr

    def render_expr(self, expr):
        node = brian_ast(expr, self.variables)
        node = self.arithmetic_simplifier.render_node(node)
        node = self.render_node(node)
        return self.node_renderer.render_node(node)

    def render_node(self, node):
        '''
        Assumes that the node has already been fully processed by BrianASTRenderer
        '''
        # can we pull this out?
        if node.scalar and node.complexity>0:
            expr = self.node_renderer.render_node(node)
            if expr in self.loop_invariants:
                name = self.loop_invariants[expr]
            else:
                self.n += 1
                name = '_lio_'+str(self.n)
                self.loop_invariants[expr] = name
                self.loop_invariant_dtypes[name] = node.dtype
            # None is the expression context, we don't use it so we just set to None
            newnode = ast.Name(name, None)
            newnode.scalar = True
            newnode.dtype = node.dtype
            newnode.complexity = 0
            return newnode
        # otherwise, render node as usual
        return super(Simplifier, self).render_node(node)
