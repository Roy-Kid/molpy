"""Unit tests for SMARTS matcher.

Tests SmartsMatcher, Candidate, and ScoringPolicy classes.
"""

import pytest
from dataclasses import dataclass
from molpy.core.atomistic import Atomistic
from molpy.typifier.matcher import SmartsMatcher, Candidate, ScoringPolicy
from molpy.typifier.builder import PatternBuilder, quick_pattern
from molpy.typifier.predicates import is_element, bond_order, is_aromatic
from molpy.typifier.adapter import build_mol_graph


class TestCandidate:
    """Test Candidate dataclass."""
    
    def test_creation(self):
        """Test creating candidate."""
        c = Candidate(
            atom_id=123,
            atomtype="carbonyl_carbon",
            source="test_pattern",
            priority=10,
            score=5,
            pattern_size=(3, 2),
            definition_order=1
        )
        
        assert c.atom_id == 123
        assert c.atomtype == "carbonyl_carbon"
        assert c.source == "test_pattern"
        assert c.priority == 10
        assert c.score == 5
        assert c.pattern_size == (3, 2)
        assert c.definition_order == 1
    
    def test_comparison_by_priority(self):
        """Test candidates compare by priority first."""
        c1 = Candidate(1, "a", "src", 10, 5, (2, 1), 1)
        c2 = Candidate(1, "b", "src", 20, 3, (2, 1), 1)
        
        # Higher priority wins
        assert c2 < c1  # c2 has higher priority, sorts earlier
    
    def test_comparison_by_score(self):
        """Test candidates compare by score when priority equal."""
        c1 = Candidate(1, "a", "src", 10, 5, (2, 1), 1)
        c2 = Candidate(1, "b", "src", 10, 8, (2, 1), 1)
        
        # Higher score wins
        assert c2 < c1  # c2 has higher score, sorts earlier
    
    def test_comparison_by_size(self):
        """Test candidates compare by size when priority and score equal."""
        c1 = Candidate(1, "a", "src", 10, 5, (2, 1), 1)
        c2 = Candidate(1, "b", "src", 10, 5, (3, 2), 1)
        
        # Larger size wins
        assert c2 < c1  # c2 has larger size, sorts earlier
    
    def test_comparison_by_definition_order(self):
        """Test candidates compare by definition order when all else equal."""
        c1 = Candidate(1, "a", "src", 10, 5, (2, 1), 1)
        c2 = Candidate(1, "b", "src", 10, 5, (2, 1), 2)
        
        # Later definition wins (higher number)
        assert c2 < c1  # c2 has later definition, sorts earlier
    
    def test_equality(self):
        """Test candidate equality."""
        c1 = Candidate(1, "a", "src", 10, 5, (2, 1), 1)
        c2 = Candidate(1, "a", "src", 10, 5, (2, 1), 1)
        
        assert c1 == c2
    
    def test_total_ordering(self):
        """Test candidates have total ordering."""
        candidates = [
            Candidate(1, "a", "src", 5, 3, (2, 1), 3),
            Candidate(1, "b", "src", 10, 2, (1, 0), 1),
            Candidate(1, "c", "src", 10, 5, (3, 2), 2),
            Candidate(1, "d", "src", 10, 5, (3, 2), 1),
        ]
        
        sorted_candidates = sorted(candidates)
        
        # Expected order after sort (best to worst):
        # c: priority=10, score=5, size=5, order=2 (later definition wins)
        # d: priority=10, score=5, size=5, order=1
        # b: priority=10, score=2, size=1, order=1
        # a: priority=5
        assert sorted_candidates[0].atomtype == "c"
        assert sorted_candidates[1].atomtype == "d"
        assert sorted_candidates[2].atomtype == "b"
        assert sorted_candidates[3].atomtype == "a"


class TestScoringPolicy:
    """Test ScoringPolicy class."""
    
    def test_default_method(self):
        """Test default scoring method."""
        # Create a simple pattern
        pattern = quick_pattern("test", "C")
        
        score = ScoringPolicy.default(pattern)
        
        # Should return non-negative score
        assert score >= 0
    
    def test_custom_method(self):
        """Test custom scoring method with weights."""
        pattern = quick_pattern("test", "C")
        
        score1 = ScoringPolicy.custom(pattern, vertex_weight=1.0, edge_weight=1.5)
        score2 = ScoringPolicy.custom(pattern, vertex_weight=2.0, edge_weight=3.0)
        
        # Higher weights should give higher scores
        assert score2 >= score1
    
    def test_scoring_specificity(self):
        """Test scoring reflects pattern specificity."""
        # Simple pattern: just element
        p1 = quick_pattern("simple", "C")
        
        # More specific pattern: element + charge
        from molpy.typifier.predicates import charge
        builder = PatternBuilder("specific")
        builder.add_vertex([is_element("C"), charge(0)])
        p2 = builder.build()
        
        score1 = ScoringPolicy.default(p1)
        score2 = ScoringPolicy.default(p2)
        
        # More specific pattern should have higher score
        assert score2 > score1


class TestSmartsMatcher:
    """Test SmartsMatcher class."""
    
    def test_initialization(self):
        """Test matcher initialization."""
        patterns = [quick_pattern("test", "C")]
        matcher = SmartsMatcher(patterns)
        
        assert len(matcher.patterns) == 1
        assert matcher.scoring is not None
    
    def test_custom_scoring(self):
        """Test matcher with custom scoring policy."""
        patterns = [quick_pattern("test", "C")]
        scoring = ScoringPolicy()
        matcher = SmartsMatcher(patterns, scoring=scoring)
        
        assert matcher.scoring is scoring
    
    def test_find_candidates_single_pattern(self):
        """Test finding candidates with single pattern."""
        mol = Atomistic()
        c = mol.add_atom(element="C")
        
        patterns = [quick_pattern("carbon", "C")]
        matcher = SmartsMatcher(patterns)
        
        graph, vs_to_atomid, _ = build_mol_graph(mol)
        candidates = matcher.find_candidates(graph, vs_to_atomid)
        
        assert len(candidates) == 1
        assert candidates[0].atomtype == "carbon"
    
    def test_find_candidates_no_matches(self):
        """Test finding candidates with no matches."""
        mol = Atomistic()
        mol.add_atom(element="C")
        
        patterns = [quick_pattern("oxygen", "O")]
        matcher = SmartsMatcher(patterns)
        
        graph, vs_to_atomid, _ = build_mol_graph(mol)
        candidates = matcher.find_candidates(graph, vs_to_atomid)
        
        assert len(candidates) == 0
    
    def test_find_candidates_multiple_matches(self):
        """Test finding candidates with multiple matches."""
        mol = Atomistic()
        c1 = mol.add_atom(element="C")
        c2 = mol.add_atom(element="C")
        
        patterns = [quick_pattern("carbon", "C")]
        matcher = SmartsMatcher(patterns)
        
        graph, vs_to_atomid, _ = build_mol_graph(mol)
        candidates = matcher.find_candidates(graph, vs_to_atomid)
        
        assert len(candidates) == 2
    
    def test_resolve_no_conflicts(self):
        """Test resolve with no conflicts."""
        mol = Atomistic()
        c = mol.add_atom(element="C")
        
        patterns = [quick_pattern("carbon", "C")]
        matcher = SmartsMatcher(patterns)
        
        graph, vs_to_atomid, _ = build_mol_graph(mol)
        candidates = matcher.find_candidates(graph, vs_to_atomid)
        typing = matcher.resolve(candidates)
        
        assert len(typing) == 1
        assert typing[id(c)] == "carbon"
    
    def test_resolve_with_conflicts_priority(self):
        """Test resolve uses priority to break ties."""
        mol = Atomistic()
        c = mol.add_atom(element="C", charge=0)
        
        # Lower priority pattern
        p1 = quick_pattern("any_carbon", "C", priority=5)
        
        # Higher priority pattern
        p2 = quick_pattern("neutral_carbon", "C", priority=10, charge=0)
        
        matcher = SmartsMatcher([p1, p2])
        graph, vs_to_atomid, _ = build_mol_graph(mol)
        candidates = matcher.find_candidates(graph, vs_to_atomid)
        typing = matcher.resolve(candidates)
        
        # Higher priority should win
        assert typing[id(c)] == "neutral_carbon"
    
    def test_resolve_with_conflicts_specificity(self):
        """Test resolve uses specificity when priority equal."""
        mol = Atomistic()
        c = mol.add_atom(element="C", charge=0)
        
        # Less specific
        p1 = quick_pattern("carbon", "C")
        
        # More specific (element + charge)
        p2 = quick_pattern("neutral_carbon", "C", charge=0)
        
        matcher = SmartsMatcher([p1, p2])
        graph, vs_to_atomid, _ = build_mol_graph(mol)
        candidates = matcher.find_candidates(graph, vs_to_atomid)
        typing = matcher.resolve(candidates)
        
        # More specific should win
        assert typing[id(c)] == "neutral_carbon"
    
    def test_resolve_with_conflicts_size(self):
        """Test resolve uses size when priority and specificity equal."""
        mol = Atomistic()
        c1 = mol.add_atom(element="C")
        c2 = mol.add_atom(element="C")
        mol.add_bond(c1, c2, order=1)
        
        # Single atom
        p1 = quick_pattern("carbon", "C")
        
        # Two atoms (larger)
        builder = PatternBuilder("c_c")
        v1 = builder.add_vertex([is_element("C")])
        v2 = builder.add_vertex([is_element("C")])
        builder.add_edge(v1, v2, [bond_order(1)])
        builder.set_target_vertices([v1])
        p2 = builder.build()
        
        matcher = SmartsMatcher([p1, p2])
        graph, vs_to_atomid, _ = build_mol_graph(mol)
        candidates = matcher.find_candidates(graph, vs_to_atomid)
        typing = matcher.resolve(candidates)
        
        # Larger pattern should win for c1
        assert typing[id(c1)] == "c_c"
    
    def test_resolve_deterministic_order(self):
        """Test resolve is deterministic."""
        mol = Atomistic()
        c = mol.add_atom(element="C")
        
        p1 = quick_pattern("p1", "C")
        p2 = quick_pattern("p2", "C")
        
        matcher = SmartsMatcher([p1, p2])
        graph, vs_to_atomid, _ = build_mol_graph(mol)
        candidates = matcher.find_candidates(graph, vs_to_atomid)
        
        typing1 = matcher.resolve(candidates)
        typing2 = matcher.resolve(candidates)
        
        # Should be same result
        assert typing1 == typing2
    
    def test_explain_mode(self):
        """Test explain mode returns conflict info."""
        mol = Atomistic()
        c = mol.add_atom(element="C")
        
        p1 = quick_pattern("p1", "C")
        p2 = quick_pattern("p2", "C")
        
        matcher = SmartsMatcher([p1, p2])
        graph, vs_to_atomid, _ = build_mol_graph(mol)
        candidates = matcher.find_candidates(graph, vs_to_atomid)
        
        explain = matcher.explain(candidates)
        
        assert explain is not None
        assert id(c) in explain
        assert len(explain[id(c)]['candidates']) == 2
    
    def test_match_single_atom(self):
        """Test matching single atom pattern."""
        mol = Atomistic()
        c = mol.add_atom(element="C", charge=0)
        
        patterns = [quick_pattern("carbon", "C")]
        matcher = SmartsMatcher(patterns)
        
        graph, vs_to_atomid, _ = build_mol_graph(mol)
        candidates = matcher.find_candidates(graph, vs_to_atomid)
        typing = matcher.resolve(candidates)
        
        assert typing[id(c)] == "carbon"
    
    def test_match_with_charge(self):
        """Test matching with charge constraint."""
        mol = Atomistic()
        c_neutral = mol.add_atom(element="C", charge=0)
        c_positive = mol.add_atom(element="C", charge=1)
        
        patterns = [quick_pattern("c_plus", "C", charge=1)]
        matcher = SmartsMatcher(patterns)
        
        graph, vs_to_atomid, _ = build_mol_graph(mol)
        candidates = matcher.find_candidates(graph, vs_to_atomid)
        typing = matcher.resolve(candidates)
        
        assert id(c_neutral) not in typing
        assert typing[id(c_positive)] == "c_plus"
    
    def test_match_bond_order(self):
        """Test matching with bond order constraint."""
        mol = Atomistic()
        c1 = mol.add_atom(element="C")
        c2 = mol.add_atom(element="C")
        o = mol.add_atom(element="O")
        
        mol.add_bond(c1, c2, order=1)
        mol.add_bond(c1, o, order=2)
        
        # Match C with double bond to O
        builder = PatternBuilder("carbonyl")
        c_v = builder.add_vertex([is_element("C")])
        o_v = builder.add_vertex([is_element("O")])
        builder.add_edge(c_v, o_v, [bond_order(2)])
        builder.set_target_vertices([c_v])
        pattern = builder.build()
        
        matcher = SmartsMatcher([pattern])
        graph, vs_to_atomid, _ = build_mol_graph(mol)
        candidates = matcher.find_candidates(graph, vs_to_atomid)
        typing = matcher.resolve(candidates)
        
        assert typing[id(c1)] == "carbonyl"
        assert id(c2) not in typing
    
    def test_match_aromatic(self):
        """Test matching aromatic atoms."""
        mol = Atomistic()
        c_aromatic = mol.add_atom(element="C", is_aromatic=True)
        c_aliphatic = mol.add_atom(element="C", is_aromatic=False)
        
        patterns = [quick_pattern("aromatic_c", "C", is_aromatic=True)]
        matcher = SmartsMatcher(patterns)
        
        graph, vs_to_atomid, _ = build_mol_graph(mol)
        candidates = matcher.find_candidates(graph, vs_to_atomid)
        typing = matcher.resolve(candidates)
        
        assert typing[id(c_aromatic)] == "aromatic_c"
        assert id(c_aliphatic) not in typing
    
    def test_match_ring(self):
        """Test matching ring atoms."""
        mol = Atomistic()
        ring_atoms = [mol.add_atom(element="C") for _ in range(6)]
        chain_atom = mol.add_atom(element="C")
        
        # Create ring
        for i in range(6):
            j = (i + 1) % 6
            mol.add_bond(ring_atoms[i], ring_atoms[j], order=1)
        
        patterns = [quick_pattern("ring_c", "C", in_ring=True)]
        matcher = SmartsMatcher(patterns)
        
        graph, vs_to_atomid, _ = build_mol_graph(mol)
        candidates = matcher.find_candidates(graph, vs_to_atomid)
        typing = matcher.resolve(candidates)
        
        # All ring atoms should match
        for atom in ring_atoms:
            assert typing[id(atom)] == "ring_c"
        
        # Chain atom should not match
        assert id(chain_atom) not in typing
    
    def test_match_hybridization(self):
        """Test matching with hybridization constraint."""
        mol = Atomistic()
        sp2_c = mol.add_atom(element="C", hyb=2)
        sp3_c = mol.add_atom(element="C", hyb=3)
        
        patterns = [quick_pattern("sp2_c", "C", hyb=2)]
        matcher = SmartsMatcher(patterns)
        
        graph, vs_to_atomid, _ = build_mol_graph(mol)
        candidates = matcher.find_candidates(graph, vs_to_atomid)
        typing = matcher.resolve(candidates)
        
        assert typing[id(sp2_c)] == "sp2_c"
        assert id(sp3_c) not in typing
    
    def test_performance_large_molecule(self):
        """Test performance on larger molecule."""
        import time
        
        # Create a molecule with ~200 atoms
        mol = Atomistic()
        atoms = []
        for i in range(200):
            atoms.append(mol.add_atom(element="C" if i % 2 == 0 else "O"))
        
        # Add some bonds
        for i in range(199):
            mol.add_bond(atoms[i], atoms[i+1], order=1)
        
        p_c = quick_pattern("carbon", "C")
        p_o = quick_pattern("oxygen", "O")
        matcher = SmartsMatcher([p_c, p_o])
        
        graph, vs_to_atomid, _ = build_mol_graph(mol)
        
        start = time.time()
        candidates = matcher.find_candidates(graph, vs_to_atomid)
        typing = matcher.resolve(candidates)
        elapsed = time.time() - start
        
        assert elapsed < 0.2  # Should complete in < 200ms
        assert len(typing) == 200


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
