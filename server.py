import web
from web import form
from clearer import calculate_burden
import math
import json
import os

render = web.template.render(os.path.dirname(os.path.abspath(__file__)) + '/templates/', base='layout')

#
# this is the 'URL router'
# it says which URLS go to which views
#
urls = (
    '/', 'index', 

    '/about', 'about', 
    '/contact', 'contact', 
    '/error', 'error',

    # would uncommment this next line to add a login view
    # '/login', 'login',
)

# this line creates the application object
# just no better place to put it, in a big app we'd 
# have better code structure
app = web.application(urls, globals(), autoreload=False)

#
# This is the form object that is renderdd on the webpage
# Here is an example of the form we used as a guide: http://webpy.org/form#example
#
form_comparators=['All CNS','Psychiatry','Cognition','Other Neurologic','Fall Risk','Custom']

myform = form.Form(
    form.Textbox('Index_drug',form.Validator('Field cannot be empty', lambda x: not x=="hello")),
    form.Dropdown('Comparator',form_comparators,value=0),
    form.Textarea('Druglist',value='None'),
    form.Checkbox('Option_1',value="0"),
    form.Checkbox('Option_2',value="0"),
  ) 

# this is one "view"
# it has GET and POST methods, which do different things
# GET is for just viewing the form, and POST is for submitting data
class index: 
    def GET(self): 

        # on a get request, 1) create the form object 
        # 2) render the page, passing the form object to the template
        # (it infers the basicform.html from the name of the basicform method) 
        form = myform()
        return render.basicform(form)

    def POST(self): 

        # post requests are more complicated, they have 2 possible paths...
        # either the form is valid or not, reember that the view is different
        form = myform() 

        # if it's not valid, do the same as the get request, just return the form
        # however, in this case, the data user submitted is *still bound to the form*, so 
        # results will pre-populate 
        # the form.render() method includes error text AND fills out previous results
        if not form.validates(): 
            return render.basicform(form)
        
        # if form IS valid, let's render a new template that displays the reuslt
        # remember that in the demo they only returned a string
        else:
        
            variables = {
            
                'Index_drug': str(form['Index_drug'].value).upper(),
                'Comparator': str(form['Comparator'].value),
                'Druglist': str(form['Druglist'].value).upper(),           
                'Option_1': (1 if form['Option_1'].checked else 0),
                'Option_2': (1 if form['Option_2'].checked else 0),
            
            }
            
            
            result=calculate_burden.check_medlist(variables)

            if result['matched_drugs']==[] : return render.error(variables['Index_drug'],variables['Comparator'])
            else: return render.results(result['matched_drugs'],result['list_by_drug'],result['list_by_ae'],result['annotation_by_drug'],result['ae_score'],result['drug_score'],result['ae_total'],result['mods_p450'],result['subs_p450'])

class error:
    def GET(self):
        return render.error(variables['Index_drug'],variables['Comparator'])
        
class about: 
    def GET(self): 
        return render.about()

class contact: 
    def GET(self): 
        return render.contact()

if __name__=="__main__":
    web.internalerror = web.debugerror
    app.run()
application = app.wsgifunc()
