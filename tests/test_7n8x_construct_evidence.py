import pytest
from scripts.collect_7n8x_construct_evidence import DOI, paragraphs


def test_publication_identification_and_focused_reconstruction_evidence_are_preserved():
    raw = f'''<article><front><article-id pub-id-type="doi">{DOI}</article-id></front>
    <body><p>We study full-length TcdB.</p><p>Local refinement of the central region.</p></body>
    <back><p>7N8X and 7N9Y identify the focused and full models.</p></back></article>'''.encode()
    assert len(paragraphs(raw)) == 3
    assert "7N8X" in paragraphs(raw)[-1]


def test_other_article_or_no_explicit_entry_binding_is_not_adopted():
    with pytest.raises(ValueError, match="DOI differs"):
        paragraphs(b'<article><article-id pub-id-type="doi">other</article-id><p>7N8X</p></article>')
    with pytest.raises(ValueError, match="explicitly identify"):
        paragraphs(f'<article><article-id pub-id-type="doi">{DOI}</article-id><p>central region</p></article>'.encode())
