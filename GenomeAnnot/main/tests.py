from django.test import TestCase, Client
from django.urls import reverse


class ViewTest(TestCase):
    def test_main_annotate_view(self):
        # Create a client
        client = Client()

        # get request to the annotate view
        url = reverse("main:annotate")
        response = client.get(url)

        # Check 200 status code
        self.assertEqual(response.status_code, 200)

        # Check template
        self.assertTemplateUsed(response, "main/main_annotate.html")

        # Check the value of active_tab
        self.assertEqual(response.context["active_tab"], "annotate")

    def test_main_blast_view(self):
        # Create a client
        client = Client()

        # get request to the blast view
        url = reverse("main:blast")
        response = client.get(url)

        # Check 200 status code
        self.assertEqual(response.status_code, 200)

        # Check template
        self.assertTemplateUsed(response, "main/main_blast.html")

        # Check the value of active_tab
        self.assertEqual(response.context["active_tab"], "blast")

    def test_main_explore_view(self):
        # Create a client
        client = Client()

        # get request to the explore view
        url = reverse("main:explore")
        response = client.get(url)

        # Check 200 status code
        self.assertEqual(response.status_code, 200)

        # Check template
        self.assertTemplateUsed(response, "main/main_explore.html")

        # Check the value of active_tab
        self.assertEqual(response.context["active_tab"], "explore")

    def test_main_validate_view(self):
        # Create a client
        client = Client()

        # get request to the validate view
        url = reverse("main:validate")
        response = client.get(url)

        # Check 200 status code
        self.assertEqual(response.status_code, 200)

        # Check template
        self.assertTemplateUsed(response, "main/main_validate.html")

        # Check the value of active_tab
        self.assertEqual(response.context["active_tab"], "validate")
