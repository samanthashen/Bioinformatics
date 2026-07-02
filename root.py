# -*- coding: utf-8 -*-
"""Main Controller"""

from tg import expose, flash, require, url, lurl
from tg import request, redirect, tmpl_context
from tg.i18n import ugettext as _, lazy_ugettext as l_
from tg.exceptions import HTTPFound

from digest.lib.base import BaseController
from digest.controllers.error import ErrorController

__all__ = ['RootController']

import tw2.forms as twf

from tg import validate
from formencode import validators, compound, schema


class SearchFormValidator(schema.Schema):
    query = compound.All(validators.Regex(r'^ *[GALMFWKQESPVICYHRNDTgalmfwkqespvicyhrndt]* *$'),
                         validators.PlainText(),
                         validators.String(min=3,strip=True))
    mode  = validators.OneOf(["Lys C", "Lys N", "CNBr", "Arg C",
                        "Asp N", "Trypsin",
                        "Trypsin (higher specificity)"])
    maxmissedcleave = validators.Int(min=0)
    minlen = validators.Int(min=0)
    maxlen = validators.Int(min=1)
    minwt = validators.Int(min=0)
    maxwt = validators.Int(min=1)

class SearchForm(twf.Form):
    class child(twf.TableLayout):
        query = twf.TextField(label="Search Term")
        mode = twf.SingleSelectField(label="Protein to Digest",
                                     options=["Lys C", "Lys N", "CNBr", "Arg C",
                                              "Asp N", "Trypsin", "Trypsin (even before P)",
                                              "Trypsin (higher specificity)"],
                                     prompt_text=None)
        maxmissedcleave = twf.TextField(label="Maximum number of missed cleavages")
        minlen = twf.TextField(label="Min. Length")
        maxlen = twf.TextField(label="Max. Length")
        minwt = twf.TextField(label="Min. Weight")
        maxwt = twf.TextField(label="Max. Weight")
        css_class = 'table'
        attrs = {'style': 'width: 600px;'}
    action = '/digest'
    submit = twf.SubmitButton(value="Search")
    validator = SearchFormValidator


from proj import proj

class RootController(BaseController):

    @expose('digest.templates.index')
    def index(self, **kwargs):
        """Handle the front-page."""
        return dict(title = "Protein Digest",
                    form= SearchForm)

    @expose('digest.templates.digest')
    @validate(SearchForm, error_handler=index)
    def digest(self, query, mode, maxmissedcleave, minlen,
                maxlen, minwt, maxwt):
        print(self, query, mode, maxmissedcleave, minlen,
                maxlen, minwt, maxwt)
        t = proj(query, mode, maxmissedcleave, minlen,
                 maxlen, minwt, maxwt)
        print(t) 
        title = "Protein Digest"
        return dict(title = title,  digest=t)

