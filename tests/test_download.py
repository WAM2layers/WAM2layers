import pytest

from wam2layers.download import get_4tu_catalog, get_4tu_id, parse_doi


@pytest.mark.parametrize(
    "doi,prefix,suffix",
    [
        ("https://doi.org/10.5281/zenodo.11093876", "10.5281", "zenodo.11093876"),
        ("http://doi.org/10.5281/zenodo.11093876", "10.5281", "zenodo.11093876"),
        ("doi.org/10.5281/zenodo.11093876", "10.5281", "zenodo.11093876"),
        ("10.5281/zenodo.11093876", "10.5281", "zenodo.11093876"),
        (
            "https://doi.org/10.4121/f9572240-f179-4338-9e1b-82c5598529e2",
            "10.4121",
            "f9572240-f179-4338-9e1b-82c5598529e2",
        ),
        (
            "http://doi.org/10.4121/f9572240-f179-4338-9e1b-82c5598529e2",
            "10.4121",
            "f9572240-f179-4338-9e1b-82c5598529e2",
        ),
        (
            "doi.org/10.4121/f9572240-f179-4338-9e1b-82c5598529e2",
            "10.4121",
            "f9572240-f179-4338-9e1b-82c5598529e2",
        ),
        (
            "10.4121/f9572240-f179-4338-9e1b-82c5598529e2",
            "10.4121",
            "f9572240-f179-4338-9e1b-82c5598529e2",
        ),
    ],
)
def test_parse_doi(doi, prefix, suffix):
    pref, suff = parse_doi(doi)
    assert pref == prefix
    assert suff == suffix


@pytest.mark.parametrize(
    "doi,id_4tu",
    [
        (
            "https://doi.org/10.4121/bbe10a2a-39dc-4098-a69f-0f677d06ecdd.v1",
            "bbe10a2a-39dc-4098-a69f-0f677d06ecdd",
        ),
        (
            "doi.org/10.4121/bbe10a2a-39dc-4098-a69f-0f677d06ecdd.v1",
            "bbe10a2a-39dc-4098-a69f-0f677d06ecdd",
        ),
        (
            "10.4121/bbe10a2a-39dc-4098-a69f-0f677d06ecdd.v1",
            "bbe10a2a-39dc-4098-a69f-0f677d06ecdd",
        ),
        (
            "https://doi.org/10.4121/bbe10a2a-39dc-4098-a69f-0f677d06ecdd",
            "bbe10a2a-39dc-4098-a69f-0f677d06ecdd",
        ),
    ],
)
def test_get_4tu_id(doi, id_4tu):
    result = get_4tu_id(doi)
    assert result == id_4tu


def test_get_4tu_catalog():
    test_id = "bbe10a2a-39dc-4098-a69f-0f677d06ecdd"
    result = get_4tu_catalog(test_id)
    assert result is not None
    assert result.endswith(".xml")
