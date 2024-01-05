from django.test import TestCase, Client
from django.urls import reverse


class ExploreViewTest(TestCase):
    def test_main_validate_view(self):
        # Create a client
        client = Client()

        # get request to the validate view
        url = reverse("validate:index")
        response = client.get(url)

        # Check 200 status code
        self.assertEqual(response.status_code, 200)

        # Check template
        self.assertTemplateUsed(response, "validate/main_validate.html")

        # Check the value of active_tab
        self.assertEqual(response.context["active_tab"], "validate")
