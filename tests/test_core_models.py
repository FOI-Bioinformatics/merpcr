#!/usr/bin/env python3
"""
Tests for core data models.
"""

import unittest
import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from merpcr.core.models import STSRecord, FASTARecord, STSHit, ThreadData


@pytest.mark.unit
class TestSTSRecord(unittest.TestCase):
    """Tests for STSRecord data model."""
    
    def test_sts_record_creation(self):
        """Test basic STS record creation."""
        sts = STSRecord(
            id="TEST001",
            primer1="ATCGATCG",
            primer2="CGATCGAT", 
            pcr_size=150,
            alias="Test STS",
            offset=10,
            hash_offset=2,
            direct='+',
            ambig_primer=0
        )
        
        self.assertEqual(sts.id, "TEST001")
        self.assertEqual(sts.primer1, "ATCGATCG")
        self.assertEqual(sts.primer2, "CGATCGAT")
        self.assertEqual(sts.pcr_size, 150)
        self.assertEqual(sts.alias, "Test STS")
        self.assertEqual(sts.offset, 10)
        self.assertEqual(sts.hash_offset, 2)
        self.assertEqual(sts.direct, '+')
        self.assertEqual(sts.ambig_primer, 0)
    
    def test_sts_record_defaults(self):
        """Test STS record with default values."""
        sts = STSRecord(
            id="TEST002",
            primer1="GGCCTTAA",
            primer2="TTAAGGCC",
            pcr_size=200
        )
        
        self.assertEqual(sts.alias, "")
        self.assertEqual(sts.offset, 0)
        self.assertEqual(sts.hash_offset, 0)
        self.assertEqual(sts.direct, '+')
        self.assertEqual(sts.ambig_primer, 0)


@pytest.mark.unit
class TestFASTARecord(unittest.TestCase):
    """Tests for FASTARecord data model."""
    
    def test_fasta_record_with_label(self):
        """Test FASTA record with explicit label."""
        fasta = FASTARecord(
            defline=">seq1 Test sequence",
            sequence="ATCGATCGATCG",
            label="seq1"
        )
        
        self.assertEqual(fasta.defline, ">seq1 Test sequence")
        self.assertEqual(fasta.sequence, "ATCGATCGATCG")
        self.assertEqual(fasta.label, "seq1")
    
    def test_fasta_record_auto_label(self):
        """Test automatic label extraction from defline."""
        fasta = FASTARecord(
            defline=">chromosome1 Homo sapiens chromosome 1",
            sequence="AAAATTTTCCCCGGGG"
        )
        
        self.assertEqual(fasta.label, "chromosome1")
    
    def test_fasta_record_no_arrow(self):
        """Test label extraction without > character."""
        fasta = FASTARecord(
            defline="seq_abc Description here",
            sequence="ATCG"
        )
        
        self.assertEqual(fasta.label, "seq_abc")
    
    def test_fasta_record_complex_defline(self):
        """Test complex defline parsing."""
        fasta = FASTARecord(
            defline=">gi|123456|gb|ABC123.1| Homo sapiens test sequence",
            sequence="GATTACA"
        )
        
        self.assertEqual(fasta.label, "gi|123456|gb|ABC123.1|")


@pytest.mark.unit
class TestSTSHit(unittest.TestCase):
    """Tests for STSHit data model."""
    
    def test_sts_hit_creation(self):
        """Test STS hit creation."""
        sts = STSRecord(
            id="TEST001",
            primer1="ATCG",
            primer2="CGAT",
            pcr_size=100
        )
        
        hit = STSHit(pos1=100, pos2=200, sts=sts)
        
        self.assertEqual(hit.pos1, 100)
        self.assertEqual(hit.pos2, 200)
        self.assertEqual(hit.sts, sts)
        self.assertEqual(hit.sts.id, "TEST001")


@pytest.mark.unit
class TestThreadData(unittest.TestCase):
    """Tests for ThreadData model."""
    
    def test_thread_data_creation(self):
        """Test thread data creation."""
        thread_data = ThreadData(
            thread_id=0,
            sequence="ATCGATCG",
            offset=1000,
            length=500
        )
        
        self.assertEqual(thread_data.thread_id, 0)
        self.assertEqual(thread_data.sequence, "ATCGATCG")
        self.assertEqual(thread_data.offset, 1000)
        self.assertEqual(thread_data.length, 500)
        self.assertEqual(len(thread_data.hits), 0)
    
    def test_thread_data_with_hits(self):
        """Test thread data with hits."""
        sts = STSRecord(id="TEST", primer1="AT", primer2="GC", pcr_size=50)
        hit = STSHit(pos1=10, pos2=60, sts=sts)
        
        thread_data = ThreadData(
            thread_id=1,
            sequence="ATCGATCG",
            offset=0,
            length=8,
            hits=[hit]
        )
        
        self.assertEqual(len(thread_data.hits), 1)
        self.assertEqual(thread_data.hits[0], hit)


if __name__ == "__main__":
    unittest.main()