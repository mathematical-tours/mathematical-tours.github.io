"""Tests for mathematical edge cases that a prose-only converter can lose."""
import unittest
from convert import braced,tex_math,Preprocessor,web_id
class ConversionChecks(unittest.TestCase):
    def test_nested_tex_groups_and_escaped_braces(self):
        body,end=braced(r'{a_{k}+\{x\}+{b}}tail',0)
        self.assertEqual(body,r'a_{k}+\{x\}+{b}');self.assertEqual(r'{a_{k}+\{x\}+{b}}tail'[end:],'tail')
    def test_accents_keep_expanded_macro_as_single_argument(self):
        self.assertEqual(tex_math(r'\hat \mathbb{Z}'),r'\hat{\mathbb{Z}}')
        self.assertEqual(tex_math(r'\bar\mathbb{R}\times\tilde\mathcal{R}'),r'\bar{\mathbb{R}}\times\tilde{\mathcal{R}}')
    def test_named_equations_do_not_advance_the_numeric_counter(self):
        pre=Preprocessor('fixture',4,[])
        out=pre.prepare(r'\eql{a=b}\begin{align}\tag{$C_1$}h&=1\\c&=d\end{align}\eql{x=y}')
        self.assertEqual(pre.equation,3);self.assertIn(r'\tag{$C_1$}',out);self.assertIn(r'\tag{4.3}',out)
    def test_apostrophe_label_has_stable_web_target(self):
        self.assertEqual(web_id("eq-proox-x'rad-nonstrong-1"),'eq-proox-x-rad-nonstrong-1')
if __name__=='__main__':unittest.main()
